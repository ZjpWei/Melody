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
  library('ConQuR')
  library('ALDEx2')
  library("ANCOMBC")
  
  # General ----
  cat('Starting model building script\n')
  start.time <- proc.time()[1]
  
  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("The analysis tag needs to be provided! Exiting...\n")
  }
  # arg1: replicate number:
  # from 1 to 100
  s <- as.numeric(args[1])
  
  # arg2: simulation scenario:
  # signature sparsity: "Sig_number"
  # signature effect direction: "Sig_effdir"
  # signature effect size: "Sig_effsz"
  # case/control sequence depth unevenness: "Sig_depth"
  scenario <- args[2]
  
  # arg3: factor for this scenario:
  loc <- args[3]
  
  # parameters
  data.loc <- paste0("./Simulation/", scenario, "/", loc, "/")
  tag <- "Models_ConQuR_batch_corrected"
  if(!dir.exists(paste0(data.loc, tag))){
    dir.create(paste0(data.loc, tag))
  }
  # Load data
  load(paste0(data.loc, "prepare_data/data.rel.", as.character(s), ".Rdata"))
  L <- length(data.rel)
  K <- ncol(data.rel[[1]]$Y)
  
  count <- NULL
  X.pool <- NULL
  Study <- NULL
  Group <- NULL
  for(l in 1:L){
    Group <- c(Group, as.character(data.rel[[l]]$X))
    Study <- c(Study, rep(as.character(l), length(data.rel[[l]]$X)))
    count <- rbind(count, data.rel[[l]]$Y)
    X.pool <- c(X.pool, data.rel[[l]]$X)
  }
  meta <- data.frame(Study = Study, Group = Group)
  meta$Study <- factor(meta$Study, levels = c("1", "2", "3", "4", "5"))
  rownames(meta) <- rownames(count)
  Y.pool <- ConQuR(tax_tab=t(rel.abd), batchid=meta$Study, covariates=meta$Group, 
                   batch_ref="1", num_core = 2)
  
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
  source("./utility/rarify.R")
  Y.rarify <- NULL
  for(l in 1:L){
    Y.rarify <- rbind(Y.rarify, .rarefy(otu.tab = data.rel[[l]]$Y, ss = 2023))
  }
  
  # Proportion normalization
  P.pool <- Y.rarify / rowSums(Y.rarify)
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
