# =============================================== #
#    Compared methods using batch-corrected data  #
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
  rm(list = ls())
  cat('Starting model building script\n')
  start.time <- proc.time()[1]

  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("The analysis tag needs to be provided! Exiting...\n")
  }
  s <- args[1]
  tag <- args[2]

  set.seed(2023)
  data.loc <- paste0("./CRC_meta_analysis/", tag, "/")
  if(!dir.exists(paste0(data.loc, "Models_batch_corrected"))){
    dir.create(paste0(data.loc, "Models_batch_corrected"))
  }
  load(paste0(data.loc, "prepare_data/data.rel.batch.", s, ".Rdata"))
  L <- length(data.rel.batch)
  K <- ncol(data.rel.batch[[1]]$Y)
  study <- NULL
  Y.pool <- NULL
  X.pool <- NULL
  for(l in 1:L){
    study <- c(study, rep(l, length(data.rel.batch[[l]]$X)) )
    Y.pool <- rbind(Y.pool, data.rel.batch[[l]]$Y)
    X.pool <- c(X.pool, data.rel.batch[[l]]$X)
  }
  ################################# CLR-LASSO ###############################
  Z.pool <- log(Y.pool + 0.5)
  clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))
  clrlasso <- cv.glmnet(x = clrx_Z, y = X.pool, nfolds = 5, family = 'binomial')
  
  # save model
  save(clrlasso, file = paste0(data.loc, "Models_batch_corrected/CLR.lasso.model.", s,".Rdata"))
  
  ################################# ALDEx2 ##################################
  # functions for ALDEx2
  source("./utility/aldex2.R")

  # combine data
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table <- data.frame(t(feature.table))
  meta.data <- data.frame(labels = factor(outcome), study = factor(study))
  colnames(feature.table) <- rownames(meta.data)

  # ALDEx2 model
  Aldex2.model <- run_aldex2(otu.tab = feature.table, meta = meta.data, formula = "labels + study")

  # save model
  save(Aldex2.model, file = paste0(data.loc, "Models_batch_corrected/Aldex2.model.", s, ".Rdata"))

  ################################# ANCOM ###################################
  # functions for ANCOM-BC
  source("./utility/ancombc.R")
  # ANCOM-BC model
  ANCOMBC.model <- ancombc.fun(feature.table = feature.table, meta = meta.data,  formula = "labels + study")

  # save model
  save(ANCOMBC.model, file = paste0(data.loc, "Models_batch_corrected/ANCOMBC.model.", s, ".Rdata"))

  ################################# BW (prop) ###############################
  source("./utility/rarify.R")
  Y.rarify <- NULL
  for(l in 1:L){
    Y.rarify <- rbind(Y.rarify, .rarefy(otu.tab = data.rel.batch[[l]]$Y, ss = 2023))
  }
  
  # Proportion normalization
  P.pool <- Y.rarify / rowSums(Y.rarify)
  pval <- NULL
  statis <- NULL
  for(k in 1:K){
    BW.data <- data.frame(y = P.pool[,k], x = factor(X.pool, levels = c("1","0")), block=as.factor(study) )
    pval <- c(pval, pvalue(wilcox_test(y ~ x | block, data = BW.data, )) )
    statis <- c(statis, statistic(wilcox_test(y ~ x | block, data = BW.data)))
  }
  # multiple testing correction
  BW.prop.model <- data.frame(p.val = pval, q.val = p.adjust(pval, method="BH"),
                              statis = statis, row.names = colnames(P.pool))

  # save model
  save(BW.prop.model, file = paste0(data.loc, "Models_batch_corrected/BW.prop.model.", s, ".Rdata"))

  
