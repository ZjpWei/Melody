  ##############################################
#                                                #
#                   Figure                       #
#                                                #
  ##############################################
  
  # Packages
  library('ggplot2')
  library('lattice')
  library('glmnet')
  
  # general 
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
  # signature effect size: "By_effsz"
  # signature sparsity: "By_taxnum"
  # signature effect direction: "By_pos_prop"
  # signature prevalence: "By_abd_prop"
  # case/control sequence depth unevenness: "By_seqdepth"
  scenario <- as.numeric(args[2])
  
  # arg3: factor for this scenario: 
  loc <- as.numeric(args[3])
  
  # Check files
  data.loc <- paste0("./Simulation/", scenario, "/", loc, "/")
  if(!dir.exists(paste0(data.loc, "AUPRC"))){
    dir.create(paste0(data.loc, "AUPRC"))
  }
  
  # Load data
  load(paste0(data.loc, "prepare_data/data.rel.1.Rdata"))
  tax.names <- colnames(data.rel[[1]]$Y)
  K <- length(tax.names)
  source("./utility/roc.R")
  
  AUPRC <- NULL
  x.label <- NULL
  method <- NULL
  linetype <- NULL
  S <- NULL
  for(s in 1:100){
    load(paste0(data.loc, "Melody_abess_wald/Melody.model.",as.character(s),".Rdata"))
    load(paste0(data.loc, "signal/signal.",as.character(s),".Rdata"))
    
    true.active <- rep(0, K)
    names(true.active) <- tax.names
    true.active[signal.idx] <- 1
    Paths <- list()
    for(ll in 1:ncol(Melody.model$coef)){
      tmp.path <- rep(0, K)
      names(tmp.path) <- tax.names
      tmp.path[rownames(Melody.model$coef)[Melody.model$coef[,ll] != 0]] <- TRUE
      Paths[[ll]] <- tmp.path
    }
    hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
    AUPRC <- c(AUPRC, hg.pr$AUC)
    S <- c(S, s)
    linetype <- c(linetype, "Original")
    method <- c(method, "Melody")
  }
  
  for(s in 1:100){
    load(paste0(data.loc, "Simulation/ANCOMBC.model.",as.character(s),".Rdata"))
    load(paste0(data.loc, "signal/signal.",as.character(s),".Rdata"))
    rank.fdr <- rank(ANCOMBC.model$res$p_val$labels1)
    true.active <- rep(0, K)
    names(true.active) <- tax.names
    true.active[signal.idx] <- 1
    Paths <- list()
    for(ll in 1:length(rank.fdr)){
      Paths[[ll]] <- rank.fdr <= ll
    }
    hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
    AUPRC <- c(AUPRC, hg.pr$AUC)
    S <- c(S, s)
    linetype <- c(linetype, "Original")
    method <- c(method, "ANCOM-BC")
  }

  for(s in 1:100){
    load(paste0(data.loc, "Simulation/Aldex2.model.",as.character(s),".Rdata"))
    load(paste0(data.loc, "signal/signal.",as.character(s),".Rdata"))
    true.active <- rep(0, K)
    names(true.active) <- tax.names
    true.active[signal.idx] <- 1
    aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval`
    rank.fdr <- rank(aldex2.fdr)
    Paths <- list()
    for(ll in 1:length(rank.fdr)){
      Paths[[ll]] <- rank.fdr <= ll
    }
    hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
    AUPRC <- c(AUPRC, hg.pr$AUC)
    S <- c(S, s)
    linetype <- c(linetype, "Original")
    method <- c(method, "ALDEx2")
  }

  for(s in 1:100){
    load(paste0(data.loc, "Simulation/BW.prop.model.",as.character(s),".Rdata"))
    load(paste0(data.loc, "signal/signal.",as.character(s),".Rdata"))
    rank.fdr <- rank(BW.prop.model$p.val)
    true.active <- rep(0, K)
    names(true.active) <- tax.names
    true.active[signal.idx] <- 1
    Paths <- list()
    for(ll in 1:length(rank.fdr)){
      Paths[[ll]] <- rank.fdr <= ll
    }
    hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
    AUPRC <- c(AUPRC, hg.pr$AUC)
    S <- c(S, s)
    linetype <- c(linetype, "Original")
    method <- c(method, "BW")
  }

  for(s in 1:100){
    load(paste0(data.loc, "Simulation/CLR.lasso.model.",as.character(s),".Rdata"))
    load(paste0(data.loc, "signal/signal.",as.character(s),".Rdata"))
    true.active <- rep(0, K)
    names(true.active) <- tax.names
    true.active[signal.idx] <- 1
    Paths <- list()
    for(ll in 1:ncol(clrlasso$glmnet.fit$beta)){
      Paths[[ll]] <- clrlasso$glmnet.fit$beta[,ll] != 0
    }
    hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
    AUPRC <- c(AUPRC, hg.pr$AUC)
    S <- c(S, s)
    linetype <- c(linetype, "Original")
    method <- c(method, "CLR-LASSO")
  }
  
  for(s in 1:100){
    load(paste0(data.loc, "Simulation_batch/ANCOMBC.model.",as.character(s),".Rdata"))
    load(paste0(data.loc, "signal/signal.",as.character(s),".Rdata"))
    rank.fdr <- rank(ANCOMBC.model$res$p_val$labels1)
    true.active <- rep(0, K)
    names(true.active) <- tax.names
    true.active[signal.idx] <- 1
    Paths <- list()
    for(ll in 1:length(rank.fdr)){
      Paths[[ll]] <- rank.fdr <= ll
    }
    hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
    AUPRC <- c(AUPRC, hg.pr$AUC)
    S <- c(S, s)
    linetype <- c(linetype, "Batch-corrected")
    method <- c(method, "ANCOM-BC")
  }
  
  for(s in 1:100){
    load(paste0(data.loc, "Simulation_batch/Aldex2.model.",as.character(s),".Rdata"))
    load(paste0(data.loc, "signal/signal.",as.character(s),".Rdata"))
    true.active <- rep(0, K)
    names(true.active) <- tax.names
    true.active[signal.idx] <- 1
    aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval`
    rank.fdr <- rank(aldex2.fdr)
    Paths <- list()
    for(ll in 1:length(rank.fdr)){
      Paths[[ll]] <- rank.fdr <= ll
    }
    hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
    AUPRC <- c(AUPRC, hg.pr$AUC)
    S <- c(S, s)
    linetype <- c(linetype, "Batch-corrected")
    method <- c(method, "ALDEx2")
  }
  
  for(s in 1:100){
    load(paste0(data.loc, "Simulation_batch/BW.prop.model.",as.character(s),".Rdata"))
    load(paste0(data.loc, "signal/signal.",as.character(s),".Rdata"))
    rank.fdr <- rank(BW.prop.model$p.val)
    true.active <- rep(0, K)
    names(true.active) <- tax.names
    true.active[signal.idx] <- 1
    Paths <- list()
    for(ll in 1:length(rank.fdr)){
      Paths[[ll]] <- rank.fdr <= ll
    }
    hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
    AUPRC <- c(AUPRC, hg.pr$AUC)
    S <- c(S, s)
    linetype <- c(linetype, "Batch-corrected")
    method <- c(method, "BW")
  }
  
  for(s in 1:100){
    load(paste0(data.loc, "Simulation_batch/CLR.lasso.model.",as.character(s),".Rdata"))
    load(paste0(data.loc, "signal/signal.",as.character(s),".Rdata"))
    true.active <- rep(0, K)
    names(true.active) <- tax.names
    true.active[signal.idx] <- 1
    Paths <- list()
    for(ll in 1:ncol(clrlasso$glmnet.fit$beta)){
      Paths[[ll]] <- clrlasso$glmnet.fit$beta[,ll] != 0
    }
    hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
    AUPRC <- c(AUPRC, hg.pr$AUC)
    S <- c(S, s)
    linetype <- c(linetype, "Batch-corrected")
    method <- c(method, "CLR-LASSO")
  }

  PRC <- data.frame(AUPRC, S = S, method = method, linetype = linetype)
  
  PRC$method <- factor(PRC$method, levels = unique(method), ordered = TRUE)
  
  save(PRC, file = paste0(data.loc, "AUPRC/PRC.Rdata"))
  
  
  