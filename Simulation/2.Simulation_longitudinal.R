# =============================================== #
#      Simulation 2. Simulation analysis on       #
#               Longitudinal data                 #
# =============================================== #
  
  # Packages ----
  library("phyloseq")
  library("coin")
  library("purrr")
  library("doParallel")
  library('tidyverse')
  library('miMeta')
  library('glmnet')
  library('ALDEx2')
  library("ANCOMBC")
  library('MMUPHin')
  
  rm(list = ls())
  # Simulation settings ----
  ## Sample size vector in simulation. 
  ## Sample size for large scenario: c(100, 120,140, 160, 180) 
  ## Sample size for small scenario: c(20, 30, 40, 50, 60) 
  Setting <- "large"
  if(Setting == "large"){
    n.sample <- c(100, 120, 140, 160, 180)
  }else if(Setting == "small"){
    n.sample <- c(20, 30, 40, 50, 60)
  }
  ## Replicate number: integer range from 1 to 100
  s <- 1
  ## Signature effect size {20, 30, 40 (default), 50, 60, 70, 80}
  Ka <- 40
  ## Signature effect direction {0.5, 0.6, 0.7(default), 0.8, 0.9, 1}
  pos.pt <- 0.7
  ## Signature effect size
  ## Sample size for large scenario: {1, 1.5, 2(default), 2.5, 3}
  ## Sample size for small scenario: {4, 4.5, 5(default), 5.5, 6}
  effect.sz <- 2
  ## Case/control sequence depth unevenness {0(default), 0.25 0.5, 0.75, 1}
  mu <-  0
  ## Directory for saving AUPRC
  data.loc <- paste0("./Simulation/", Setting, "/Longitudinal/AUPRC_Ka", Ka, "_Pos", pos.pt, "_effsz", effect.sz, "_mu", mu, "_", s,".Rdata")  
  
  # Simulate data ----
  ## random seed
  seed <- 2023
  ## Signature sparsity 
  abd.pt <- 0.5
  
  ## Data processing, remove samples with sequence depth < 2000, apply 0.2 prevalence filter.
  load("./Data/CRC_data/data/count.Rdata")
  load("./Data/CRC_data/data/meta.Rdata")
  
  meta <- as.data.frame(meta)
  study <- c("AT-CRC", "CN-CRC", "DE-CRC", "FR-CRC", "US-CRC")
  rownames(meta) <- meta$Sample_ID
  meta$Group <- factor(meta$Group, level = c("CTR", "CRC"))
  meta$Study <- factor(meta$Study, levels = study)
  sample.id.kp <- names(which(rowSums(count) >= 2000)) 
  meta <- meta[sample.id.kp,]
  count <- count[sample.id.kp,]
  pre.filter <- colMeans(count != 0) >= 0.2
  count <- count[,pre.filter]
  
  L <- length(study)
  K <- ncol(count)
  Y.study <- list()
  X.study <- list()
  for(l in 1:L){
    id.set <- meta$Study == study[l]
    condition <- meta$Group[id.set]
    Y.study[[l]] <- count[id.set,]
    X.study[[l]] <- c(condition == 'CRC') + 0
  }
  
  ## Generate signal and signature size
  source("./utility/GDM_utility.R")
  set.seed(s + seed)
  ## Randomize uneveness
  uneveness.ind <- rbernoulli(n = L, p = 0.5)
  ## Randomize signal
  ave.prop <- colMeans(count / rowSums(count))
  ## Detect most abundant taxa and less abundant taxa
  signal.most.abund <- which(ave.prop >= 1e-3)
  signal.less.abund <- which(ave.prop  < 1e-3)
  signal.m.abd <- sample(signal.most.abund)
  signal.l.abd <- sample(signal.less.abund)
  tmp.signal.idx <- c(signal.m.abd[1:round(abd.pt * Ka)], signal.l.abd[1:(Ka - round(abd.pt * Ka))])
  ## Decide signal location and direction
  signal.sig <- sample(c(rep(TRUE, round(pos.pt * Ka)), rep(FALSE, Ka - round(pos.pt*Ka))))
  signal.idx = sort(tmp.signal.idx)
  ## Add signal
  signal.add <- c()
  tmp.add <- runif(n = Ka, min = 1, max = effect.sz + 1)
  for(l in 1:L){
    signal.add <- rbind(signal.add, tmp.add)
  }
  
  ## Simulate null data
  data.null.cluster <- list()
  for(l in 1:L){
    load(paste0("./Simulation/GDM_fit/GDM.", as.character(l), ".Rdata"))
    X <- matrix(1, nrow = n.sample[l], ncol = 1)
    dmData <- r.GDM(X = X, mod = mod.gdm)
    data.null.cluster[[l]] <- list(Y = dmData, X = rep(0, nrow(dmData)))
  }
  
  # Simulate longitudinal null data
  data.null <- list()
  for(l in 1:L){
    load(paste0("./Simulation/GDM_fit/GDM.", as.character(l), ".Rdata"))
    X.base <- matrix(1, nrow = n.sample[l]/2, ncol = 1)
    cluster <- rep(1:(n.sample[l]/2), each = 2)
    dmData <- data.null.cluster[[l]]$Y
    dmData.base <- r.GDM(X = X.base, mod = mod.gdm)
    for(ll in 1:(n.sample[l]/2)){
      dmData[cluster == ll, ] <- t(t(dmData[cluster == ll, ]) + dmData.base[ll,])
    }
    Y.tmp <- Y.study[[l]]
    idx.rev <- rank(-colMeans(Y.tmp/rowSums(Y.tmp)))
    stopifnot(sum(colMeans(Y.tmp/rowSums(Y.tmp)) != 0) == ncol(dmData))
    dmData <- cbind(dmData, matrix(0, nrow = n.sample[l], ncol = K - ncol(dmData)))[,idx.rev]
    colnames(dmData) <- colnames(Y.study[[l]])
    data.null[[l]] <- list(Y = dmData / 2, X = rep(0, nrow(dmData)), cluster = cluster)
  }

  # Simulate relative abundant data
  data.rel <- list()
  for(l in 1:L){
    org.data <- Y.study[[l]]
    c.name <- colnames(org.data)
    cluster <- data.null[[l]]$cluster
    n = length(unique(cluster))
    
    # Shuffle subjects (we always assign the first half subjects to cases)
    case.idx.1 = 1:round(n/2)
    Simulate.count.1 = data.null[[l]]$Y
    Simulate.otc.1 = data.null[[l]]$X
    for(ll in 1:Ka){
      if(signal.sig[ll]){
        Simulate.count.1[cluster %in% case.idx.1, signal.idx[ll]] <- Simulate.count.1[cluster %in% case.idx.1, signal.idx[ll]] * signal.add[l,ll]
      }else{
        Simulate.count.1[!(cluster %in% case.idx.1), signal.idx[ll]]<- Simulate.count.1[!(cluster %in% case.idx.1),signal.idx[ll]] * signal.add[l,ll]
      }
    }
    Simulate.otc.1[cluster %in% case.idx.1] = 1
    Prob.abs.1 = Simulate.count.1/rowSums(Simulate.count.1)
    
    # Simulate relative abundant data
    Simulate.depth.1 <- rep(0, n.sample[l])
    if(uneveness.ind[l]){
      Simulate.depth.1[cluster %in% case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = sum(cluster %in% case.idx.1), replace = TRUE) * (mu + 1)
      Simulate.depth.1[!(cluster %in% case.idx.1)]<- sample(x = rowSums(Y.study[[l]]), size = n.sample[l] - sum(cluster %in% case.idx.1), replace = TRUE)
    }else{
      Simulate.depth.1[cluster %in% case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = sum(cluster %in% case.idx.1), replace = TRUE)
      Simulate.depth.1[!(cluster %in% case.idx.1)]<- sample(x = rowSums(Y.study[[l]]), size = n.sample[l] - sum(cluster %in% case.idx.1), replace = TRUE) * (mu + 1)
    }
    Simulate.count.rel.1 = NULL
    for(ll in 1:nrow(Simulate.count.1)){
      if(Simulate.otc.1[ll] == 1){
        Simulate.count.rel.1 <- cbind(Simulate.count.rel.1, rmultinom(1, Simulate.depth.1[ll], Prob.abs.1[ll,]))
      }else{
        Simulate.count.rel.1 <- cbind(Simulate.count.rel.1, rmultinom(1, Simulate.depth.1[ll], Prob.abs.1[ll,]))
      }
    }
    Simulate.count.rel.1 <- t(Simulate.count.rel.1)
    
    # Assign taxa and sample names
    data.rel[[l]] <- list(Y=Simulate.count.rel.1, X=Simulate.otc.1, cluster = cluster)
    sample.nm <- paste0("Cohort.",as.character(l),".s",as.character(1:n.sample[l]))
    colnames(data.rel[[l]]$Y) <- c.name
    rownames(data.rel[[l]]$Y) <- sample.nm
  }

  ## Summarize data
  signal.names <- colnames(data.rel[[1]]$Y)[signal.idx]
  
  ## Batch correction by MMUPHin
  Study <- NULL
  Group <- NULL
  rel.abd <- NULL
  for(l in 1:L){
    Group <- c(Group, as.character(data.rel[[l]]$X))
    Study <- c(Study, rep(as.character(l), length(data.rel[[l]]$X)))
    rel.abd <- rbind(rel.abd, data.rel[[l]]$Y)
  }
  meta <- data.frame(Study = Study, Group = Group)
  meta$Study <- factor(meta$Study, levels = c("1", "2", "3", "4", "5"))
  rownames(meta) <- rownames(rel.abd)
  batch_adj <- MMUPHin::adjust_batch(feature_abd = t(rel.abd),
                                     batch = "Study",
                                     covariates = "Group",
                                     data = meta)
  Y.pool <- t(batch_adj$feature_abd_adj)
  
  data.rel.mmuphin <- list()
  for(l in 1:L){
    data.rel.mmuphin[[l]] <- list(Y = Y.pool[meta$Study == as.character(l),],
                                  X = data.rel[[l]]$X,
                                  cluster = data.rel[[l]]$cluster)
  }
  
  # ============================================================================================== #
  # Analysis on original data ----
  rm(list = setdiff(ls(), c("data.rel", "data.rel.mmuphin", "signal.names", "s", "data.loc", "Ka")))
  source("./utility/roc.R")
  AUPRC <- NULL
  S <- NULL
  method <- NULL
  batch_type <- NULL
  data_type <- NULL
  
  L <- length(data.rel)  
  K <- ncol(data.rel[[1]]$Y)
  study = NULL
  Y.pool = NULL
  X.pool = NULL
  cluster <- NULL
  for(l in 1:L){
    cluster <- c(cluster, paste0(l,"_", data.rel[[l]]$cluster))
    study = c(study, rep(l, length(data.rel[[l]]$X)) )
    Y.pool = rbind(Y.pool, data.rel[[l]]$Y)
    X.pool = c(X.pool, data.rel[[l]]$X)
  }
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(study), cluster = cluster)
  colnames(feature.table) <- rownames(meta.data)

  ## ANCOM-BC2
  source("./utility/ancombc.R")
  
  ANCOMBC2.model <- ancombc.fun(feature.table = feature.table,
                                meta = meta.data,
                                formula = "labels + study",
                                adjust.method = "fdr",
                                group = NULL,
                                subject = "cluster",
                                method = "ancombc2")
  
  ## Calculate AUPRC
  rank.fdr <- rank(ANCOMBC2.model$res$p_labels1)
  true.active <- rep(0, length(ANCOMBC2.model$res$taxon))
  names(true.active) <- ANCOMBC2.model$res$taxon
  true.active[intersect(signal.names, ANCOMBC2.model$res$taxon)] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUC)
  S <- c(S, s)
  batch_type <- c(batch_type, "Original")
  data_type <- c(data_type, "Longitudinal")
  method <- c(method, "ANCOM-BC2")
  
  ## MMUPHin 
  ## function use for ANCOMBC2
  source("./utility/mmuphin.R")
  
  ## MMUPHin random effect model
  MMUPHin.model <- fit_metaAnalysis(feature.abd = feature.table,
                                    data = meta.data,
                                    test_variable = "labels",
                                    contrasts = list("1" = "1", "0" = "0"),
                                    batch_variable = "study",
                                    covariates = NULL,
                                    covariates.random = "cluster",
                                    moderator_variables = NULL)
  
  ## Calculate AUPRC
  rank.fdr <- rank(MMUPHin.model$meta_fits$pval)
  true.active <- rep(0, nrow(MMUPHin.model$meta_fits))
  names(true.active) <- MMUPHin.model$meta_fits$feature
  true.active[signal.names] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUC)
  S <- c(S, s)
  batch_type <- c(batch_type, "Original")
  data_type <- c(data_type, "Longitudinal")
  method <- c(method, "MMUPHin")

  ## Melody
  rel.abd <- list()
  covariate.interest <- list()
  cluster <- list()
  ref <- c()
  for(l in 1:length(data.rel)){
    ref <- c(ref, "Coprococcus catus [ref_mOTU_v2_4874]")
    rel.abd[[paste0("S",l)]] <- data.rel[[l]]$Y
    covariate.interest[[paste0("S",l)]] <- data.frame(disease = data.rel[[l]]$X)
    cluster[[paste0("S",l)]] <- data.rel[[l]]$cluster
  }
  names(ref) <- names(covariate.interest)

  ## Generate null model
  null.obj <- melody.null.model(rel.abd = rel.abd, ref = ref)

  ## Generate summary statistics
  summary.stat.study <- melody.get.summary(null.obj = null.obj,
                                           covariate.interest = covariate.interest,
                                           cluster = cluster)

  ## Meta-analysis
  Melody.model <- melody.meta.summary(summary.stats = summary.stat.study,
                                      tune.path = "sequence",
                                      tune.size.sequence = 1:(Ka + 100),
                                      output.best.one = FALSE)

  true.active <- rep(0, nrow(Melody.model$disease$coef))
  names(true.active) <- rownames(Melody.model$disease$coef)
  true.active[intersect(signal.names, rownames(Melody.model$disease$coef))] <- 1
  Paths <- list()
  for(ll in 1:ncol(Melody.model$disease$coef)){
    Paths[[ll]] <- Melody.model$disease$coef[,ll]!=0
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUC)
  S <- c(S, s)
  batch_type <- c(batch_type, "Original")
  data_type <- c(data_type, "Longitudinal")
  method <- c(method, "Melody")
  
  # ============================================================================================== #
  # Analysis on MMUPHin batch corrected data ----
  rm(list = setdiff(ls(), c("AUPRC", "S", "data_type", "batch_type", "method", 
                            "s", "data.rel.mmuphin", "signal.names", "data.loc")))
  
  source("./utility/roc.R")
  L <- length(data.rel.mmuphin)
  K <- ncol(data.rel.mmuphin[[1]]$Y)
  study = NULL
  Y.pool = NULL
  X.pool = NULL
  cluster <- NULL
  for(l in 1:L){
    cluster <- c(cluster, paste0(l,"_", data.rel.mmuphin[[l]]$cluster))
    study = c(study, rep(l, length(data.rel.mmuphin[[l]]$X)) )
    Y.pool = rbind(Y.pool, data.rel.mmuphin[[l]]$Y)
    X.pool = c(X.pool, data.rel.mmuphin[[l]]$X)
  }
  
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(study), cluster = cluster)
  colnames(feature.table) <- rownames(meta.data)

  # ANCOM-BC2
  source("./utility/ancombc.R")
  ANCOMBC2.model <- ancombc.fun(feature.table = feature.table,
                                meta = meta.data,
                                formula = "labels + study",
                                adjust.method = "fdr",
                                group = NULL,
                                subject = "cluster",
                                method = "ancombc2")
  
  ## Calculate AUPRC
  rank.fdr <- rank(ANCOMBC2.model$res$p_labels1)
  true.active <- rep(0, length(ANCOMBC2.model$res$taxon))
  names(true.active) <- ANCOMBC2.model$res$taxon
  true.active[intersect(signal.names, ANCOMBC2.model$res$taxon)] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUC)
  S <- c(S, s)
  batch_type <- c(batch_type, "MMUPHin")
  data_type <- c(data_type, "Longitudinal")
  method <- c(method, "ANCOM-BC2")

  ## MMUPHin
  ## function use for ANCOM
  source("./utility/mmuphin.R")

  ## MMUPHin random effect model
  MMUPHin.model <- fit_metaAnalysis(feature.abd = feature.table,
                                    data = meta.data,
                                    test_variable = "labels",
                                    contrasts = list("1" = "1", "0" = "0"),
                                    batch_variable = "study",
                                    covariates = NULL,
                                    covariates.random = "cluster",
                                    moderator_variables = NULL)
  
  # Calculate AUPRC
  rank.fdr <- rank(MMUPHin.model$meta_fits$pval)
  true.active <- rep(0, nrow(MMUPHin.model$meta_fits))
  names(true.active) <- MMUPHin.model$meta_fits$feature
  true.active[signal.names] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUC)
  S <- c(S, s)
  batch_type <- c(batch_type, "MMUPHin")
  data_type <- c(data_type, "Longitudinal")
  method <- c(method, "MMUPHin")

  # ============================================================================================== #
  # Save AUPRC ----
  PRC <- data.frame(AUPRC = AUPRC, S = S, method = method, batch_type = batch_type, data_type = data_type)
  
  save(PRC, file = data.loc)
  
