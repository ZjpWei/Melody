# =============================================== #
#     Random split: Compared methods using        #
#                 Original data                   #
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

  # Parameters
  data.loc <- "./CRC_meta_analysis/Prediction/Random_split/"
  seed <- 2023
  tag <- "Models_original"
  if(!dir.exists(paste0(data.loc, tag))){
    dir.create(paste0(data.loc, tag))
  }

  # Load data
  load(paste0(data.loc, "/prepare_data/data.rel.all.Rdata"))
  L <- length(data.rel)
  K <- ncol(data.rel[[1]]$Y)

  # sample case and control
  set.seed(s + seed)
  case.id.lst <- list()
  control.id.lst <- list()
  for(l in 1:L){
    case.id.lst[[l]] <- sample(which(data.rel[[l]]$X == 1))
    control.id.lst[[l]] <- sample(which(data.rel[[l]]$X == 0))
  }

  # Split analysis data, training data and testing data
  data.rel.analysis <- list()
  data.rel.train <- list()
  data.rel.test <- list()
  set.seed(s + seed)
  for(l in 1:L){
    case.id <- case.id.lst[[l]]
    control.id <- control.id.lst[[l]]
    case.rg <- round(c(0.6,0.8,1) * length(case.id))
    control.rg <- round(c(0.6,0.8,1) * length(control.id))
    data.rel.analysis[[l]] <- list(Y = data.rel[[l]]$Y[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]])),],
                                   X = data.rel[[l]]$X[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]]))])

    data.rel.train[[l]] <- list(Y = data.rel[[l]]$Y[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]])),],
                                X = data.rel[[l]]$X[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]]))])

    data.rel.test[[l]] <- list(Y = data.rel[[l]]$Y[sort(c(case.id[(case.rg[2]+1):case.rg[3]], control.id[(control.rg[2]+1):control.rg[3]])),],
                               X = data.rel[[l]]$X[sort(c(case.id[(case.rg[2]+1):case.rg[3]], control.id[(control.rg[2]+1):control.rg[3]]))])
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
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(Study))
  colnames(feature.table) <- rownames(meta.data)

  ## ALDEx2 will add 0.5 pseudo-count in it's own function.
  Aldex2.model <- run_aldex2(otu.tab = feature.table, meta = meta.data, formula = "labels + study")

  # save model
  save(Aldex2.model, file = paste0(data.loc, tag, "/Aldex2.model.",as.character(s),".Rdata"))

  ################################# ANCOM-BC #################################
  # function use for ANCOM-BC
  source("./utility/ancombc.R")

  # ANCOM model
  # ANCOMBC will add 1 as pseudo-count in it's own function;
  ANCOMBC.model <- ancombc.fun(feature.table = feature.table, meta = meta.data,  formula = "labels + study")

  # save model
  save(ANCOMBC.model, file = paste0(data.loc, tag, "/ANCOMBC.model.",as.character(s),".Rdata"))

  ################################# CLR-LASSO ################################
  # Add pseudo-count 0.5
  Z.pool <- log(Y.pool + 0.5)

  clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

  clrlasso <- cv.glmnet(x = clrx_Z, y = X.pool, nfolds = 5, family = 'binomial')

  save(clrlasso, file = paste0(data.loc, tag, "/CLR.lasso.model.",as.character(s),".Rdata"))

  ################################# BW (prop) ################################
  source(paste0("./utility/rarify.R"))
  Y.rarify <- NULL
  for(l in 1:L){
    Y.rarify = rbind(Y.rarify, .rarefy(otu.tab = data.rel.analysis[[l]]$Y, ss = 2023))
  }
  
  # Proportion normalization
  P.pool = Y.rarify / rowSums(Y.rarify)
  pval = NULL
  for(k in 1:K){
    BW.data = data.frame(y=P.pool[,k], x = as.factor(X.pool), block=as.factor(Study) )
    pval = c(pval, pvalue(wilcox_test(y ~ x | block, data=BW.data)) )
  }
  # multiple testing correction
  BW.prop.model <- data.frame(p.val = pval, q.val = p.adjust(pval, method="BH"),
                              row.names = colnames(P.pool))
  
  # save model
  save(BW.prop.model, file = paste0(data.loc, tag, "/BW.prop.model.",as.character(s),".Rdata"))
  