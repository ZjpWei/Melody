# =============================================== #
#           LOSO: Compared methods using          #
#               batch-corrected data              #
# =============================================== #

  # Packages ----
  library('coin')
  library('phyloseq')
  library('pracma')
  library('corpcor')
  library('glmnet')
  library('MASS')
  library('tidyverse')
  library('ALDEx2')
  library("ANCOMBC")

  # General ----
  cat('Starting model building script\n')
  start.time <- proc.time()[1]

  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("The analysis tag needs to be provided! Exiting...\n")
  }
  s <- as.numeric(args[1])
  ss <- as.numeric(args[2])

  # Parameters
  data.loc <- "./CRC_meta_analysis/Prediction/LOSO/"
  seed <- 2023
  tag <- "Models_batch_corrected"
  if(!dir.exists(paste0(data.loc, tag))){
    dir.create(paste0(data.loc, tag))
  }

  # Load data
  load(paste0(data.loc, "/prepare_data/data.rel.batch.all.Rdata"))

  # leave testing study out
  data.rel.test <- data.rel.batch[s]
  data.rel <- data.rel.batch[-s]
  L <- length(data.rel.batch)
  K <- ncol(data.rel.batch[[1]]$Y)

  # split analysis data and training data
  data.rel.analysis <- list()
  data.rel.train <- list()
  set.seed(ss + seed)
  for(l in 1:L){
    # sample case's and control's id
    case.id <- sample(which(data.rel.batch[[l]]$X == 1))
    control.id <- sample(which(data.rel.batch[[l]]$X == 0))

    case.rg <- round(c(0.75,1) * length(case.id))
    control.rg <- round(c(0.75,1) * length(control.id))
    # generate analysis data and training data
    data.rel.analysis[[l]] <- list(Y = data.rel.batch[[l]]$Y[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]])),],
                                   X = data.rel.batch[[l]]$X[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]]))])

    data.rel.train[[l]] <- list(Y = data.rel.batch[[l]]$Y[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]])),],
                                X = data.rel.batch[[l]]$X[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]]))])
  }
  study = NULL
  Y.pool = NULL
  X.pool = NULL
  for(l in 1:L){
    study = c(study, rep(l, length(data.rel.analysis[[l]]$X)) )
    Y.pool = rbind(Y.pool, data.rel.analysis[[l]]$Y)
    X.pool = c(X.pool, data.rel.analysis[[l]]$X)
  }

  ################################# ALDEx2 ###################################
  # function use for ALDEx2
  source("./utility/aldex2.R")

  # combine data
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table <- data.frame(t(feature.table))
  meta.data <- data.frame(labels = factor(outcome), study = factor(Study))
  colnames(feature.table) <- rownames(meta.data)

  ## ALDEx2 will add 0.5 pseudo-count in it's own function.
  Aldex2.model <- run_aldex2(otu.tab = feature.table, meta = meta.data, formula = "labels + study")

  # save model
  save(Aldex2.model, file = paste0(data.loc, tag, "/Aldex2.model.",as.character(s),".",as.character(ss), ".Rdata"))

  ################################# ANCOM-BC #################################
  # function use for ANCOM-BC
  source("./utility/ancombc.R")

  # ANCOM model
  # ANCOMBC will add 1 as pseudo-count in it's own function;
  ANCOMBC.model <- ancombc.fun(feature.table = feature.table, meta = meta.data,  formula = "labels + study")

  # save model
  save(ANCOMBC.model, file = paste0(data.loc, tag, "/ANCOMBC.model.",as.character(s),".",as.character(ss), ".Rdata"))

  ################################# CLR-LASSO ################################
  # Add pseudo-count 0.5
  Z.pool <- log(Y.pool + 0.5)

  clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

  clrlasso <- cv.glmnet(x = clrx_Z, y = X.pool, nfolds = 5, family = 'binomial')

  save(clrlasso, file = paste0(data.loc, tag, "/CLR.lasso.model.",as.character(s),".",as.character(ss), ".Rdata"))

  ################################# BW (prop) ################################
  source("./utility/rarify.R")
  Y.rarify <- NULL
  for(l in 1:L){
    Y.rarify = rbind(Y.rarify, .rarefy(otu.tab = data.rel.analysis[[l]]$Y, ss = 2023))
  }
  
  # Proportion normalization
  P.pool <- Y.rarify / rowSums(Y.rarify)
  pval <- NULL
  for(k in 1:K){
    BW.data <- data.frame(y=P.pool[,k], x = as.factor(X.pool), block=as.factor(Study) )
    pval <- c(pval, pvalue(wilcox_test(y ~ x | block, data=BW.data)) )
  }
  # multiple testing correction
  BW.prop.model <- data.frame(p.val = pval, q.val = p.adjust(pval, method="BH"),
                              row.names = colnames(P.pool))
  
  # save model
  save(BW.prop.model, file = paste0(data.loc, tag, "/BW.prop.model.",as.character(s),".",as.character(ss), ".Rdata"))
  