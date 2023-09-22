  ##############################################
#                                                #
#          Meta-Analysis Model Building          #
#                                                #
  ##############################################
  
  # Packages
  library('coin')
  library('phyloseq')
  library('pracma')
  library('corpcor')
  library('glmnet')
  library('MASS')
  library('tidyverse')
  library('ALDEx2')
  library("ANCOMBC")
  
  # general
  rm(list = ls())
  cat('Starting model building script\n')
  start.time <- proc.time()[1]
  
  # args = commandArgs(trailingOnly=TRUE)
  # if (length(args)==0) {
  #   stop("The analysis tag needs to be provided! Exiting...\n")
  # }
  s <- "all" #args[1]
  tag <- "CRC_all_order" #args[2]
  
  set.seed(2023)
  data.loc <- paste0("./CRC_Real/", tag, "/")
  load(paste0(data.loc, "prepare_data/data.rel.", s, ".Rdata"))
  L <- length(data.rel)
  K <- ncol(data.rel[[1]]$Y)
  study = NULL
  Y.pool = NULL
  X.pool = NULL
  for(l in 1:L){
    study = c(study, rep(l, length(data.rel[[l]]$X)) )
    Y.pool = rbind(Y.pool, data.rel[[l]]$Y)
    X.pool = c(X.pool, data.rel[[l]]$X)
  }
  # ################################# BW (prop) ###############################
  # Proportion normalization
  P.pool = Y.pool / rowSums(Y.pool)
  pval = NULL
  for(k in 1:K){
    BW.data = data.frame(y=P.pool[,k], x = as.factor(X.pool), block=as.factor(study) )
    pval = c(pval, pvalue(wilcox_test(y ~ x | block, data=BW.data)) )
  }
  # multiple testing correction
  BW.prop.model <- data.frame(p.val = pval, q.val = p.adjust(pval, method="BH"),
                              row.names = colnames(P.pool))
  
  # save model
  save(BW.prop.model, file = paste0(data.loc, "Models/BW.prop.model.", s, ".Rdata"))
  
  # ################################# ALDEx2 ##################################
  # functions for ALDEx2
  source("./utility/aldex2.R")

  # combine data
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(study))
  colnames(feature.table) <- rownames(meta.data)
  
  # ALDEx2 model
  Aldex2.model <- run_aldex2(otu.tab = feature.table, meta = meta.data, formula = "labels + study")

  # save model
  save(Aldex2.model, file = paste0(data.loc, "Models/Aldex2.model.", s, ".Rdata"))

  # ################################# ANCOM ###################################
  # functions for ANCOM-BC
  source("./utility/ancombc.R")
  # ANCOM-BC model
  ANCOMBC.model <- ancombc.fun(feature.table = feature.table, meta = meta.data,  formula = "labels + study")

  # save model
  save(ANCOMBC.model, file = paste0(data.loc, "Models/ANCOMBC.model.", s, ".Rdata"))
  
  # ################################# CLR-LASSO ###############################
  Z.pool <- log(Y.pool + 0.5)
  clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))
  clrlasso <- cv.glmnet(x = clrx_Z, y = X.pool, family = 'binomial')
  
  # save model
  save(clrlasso, file = paste0(data.loc, "Models/CLR.lasso.model.", s,".Rdata"))
