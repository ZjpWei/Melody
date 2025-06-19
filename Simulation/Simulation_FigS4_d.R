# =============================================== #
#    Simulation set (1). Simulation analysis on   #
#               Independent data                  #
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
  library('miMeta')
  library('abess')

  rm(list = ls())

  # Simulation settings ----
  n.sample <- c(100, 120, 140, 160, 180)

  ## ---- Replicate number: integer range from 1 to 100
  s <- 1

  ## ---- Signature effect size {4, 20, 40, 80, 160}
  Ka <- 4

  ## ---- Signature effect direction {0.6, 0.8, 1}
  pos.pt <- 0.6

  ## ---- Signature effect size
  ## ---- Sample size for large scenario: 2
  ## ---- Sample size for small scenario: 5
  effect.sz <- 2

  ## ---- Case/control sequence depth unevenness {0, 0.5, 1}
  mu <-  0

  ## ---- Noise level: {0.5, 1, 1.5, 2}
  B.level <-  1
  sigma_l <- 0.2 * c(1:5) ^ B.level

  ## ---- Directory for saving AUPRC
  data.loc <- paste0("./Simulation/large/Independent_FigS4d/AUPRC_Ka", Ka, "_Pos", pos.pt, "_effsz", effect.sz, "_mu", mu, "_B", B.level, "_", s,".Rdata")

  # Simulate data ----
  ## ---- random seed
  seed <- 2023
  ## ---- Signature sparsity
  abd.pt <- 0.5

  ## ---- Data processing, remove samples with sequence depth < 2000, apply 0.2 prevalence filter.
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

  ## ---- Generate signal and signature size
  source("./utility/GDM_utility.R")
  set.seed(s + seed)

  ## ---- Randomize uneveness
  uneveness.ind <- rbernoulli(n = L, p = 0.5)

  ## ---- Randomize signal
  ave.prop <- colMeans(count / rowSums(count))

  ## ---- Detect most abundant taxa and less abundant taxa
  signal.most.abund <- which(ave.prop >= 1e-3)
  signal.less.abund <- which(ave.prop  < 1e-3)
  signal.m.abd <- sample(signal.most.abund)
  signal.l.abd <- sample(signal.less.abund)
  tmp.signal.idx <- c(signal.m.abd[1:round(abd.pt * Ka)], signal.l.abd[1:(Ka - round(abd.pt * Ka))])

  ## ---- Decide signal location and direction
  signal.sig <- sample(c(rep(TRUE, round(pos.pt * Ka)), rep(FALSE, Ka - round(pos.pt*Ka))))
  signal.idx = sort(tmp.signal.idx)

  ## ---- Add signal
  signal.add <- c()
  tmp.add <- runif(n = Ka, min = 1, max = effect.sz + 1)
  for(l in 1:L){
    signal.add <- rbind(signal.add, tmp.add)
  }

  ## ---- Simulate relative abundant data
  data.rel <- list()
  for(l in 1:L){
    org.data <- Y.study[[l]]
    c.name <- colnames(org.data)
    n = n.sample[l]

    ## ---- Fit GDM with parameters
    #load(paste0("./Simulation/GDM_fit/GDM.", as.character(l), ".Rdata"))
    load(paste0("./GDM_fit/GDM.", as.character(l), ".Rdata"))
    X <- matrix(1, nrow = n, ncol = 1)
    dmData <- r.GDM(X = X, mod = mod.gdm)
    idx <- order(colMeans(org.data/rowSums(org.data)), decreasing = TRUE)
    idx.rev <- rank(-colMeans(org.data/rowSums(org.data)))

    ## ---- Stop if dimension doesn't match
    stopifnot(sum(colMeans(org.data/rowSums(org.data)) != 0) == ncol(dmData))
    dmData <- cbind(dmData, matrix(0, nrow = n.sample[l], ncol = K - ncol(dmData)))[,idx.rev]
    colnames(dmData) <- c.name

    ## ---- Simulate microbial lood
    total.load <- pmax(stats::rnbinom(n = n, mu = 1e8, size = 1), 1000)

    ## ---- Construct AA
    Simulate.count.1 <- dmData / rowSums(dmData) * total.load
    Simulate.otc.1 = rep(0, n)

    ## ---- Add differential AA signal using spike-in (we always regard the first half subjects as cases)
    case.idx.1 = 1:round(n/2)
    for(ll in 1:Ka){
      if(signal.sig[ll]){
        Simulate.count.1[case.idx.1, signal.idx[ll]] <- Simulate.count.1[case.idx.1, signal.idx[ll]] * signal.add[l,ll]
      }else{
        Simulate.count.1[-case.idx.1, signal.idx[ll]]<- Simulate.count.1[-case.idx.1,signal.idx[ll]] * signal.add[l,ll]
      }
    }
    Simulate.otc.1[case.idx.1] = 1

    ## ---- Add feature-specific processing bias.
    log.b <- rnorm(n = ncol(Simulate.count.1), mean = 0, sd = sigma_l[l])
    for(k in 1:length(log.b)){
      Simulate.count.1[,k] <- Simulate.count.1[,k] * exp(log.b[k])
    }
    Prob.abs.1 = Simulate.count.1/rowSums(Simulate.count.1)

    ## ---- Generate observed counts
    Simulate.depth.1 <- rep(0, n)
    if(uneveness.ind[l]){
      Simulate.depth.1[case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = length(case.idx.1), replace = TRUE) * (mu + 1)
      Simulate.depth.1[-case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = n - length(case.idx.1), replace = TRUE)
    }else{
      Simulate.depth.1[case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = length(case.idx.1), replace = TRUE)
      Simulate.depth.1[-case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = n - length(case.idx.1), replace = TRUE) * (mu + 1)
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

    ## ---- Assign taxa and sample names
    data.rel[[l]] <- list(Y=Simulate.count.rel.1, X=Simulate.otc.1)
    sample.nm <- paste0("Cohort.",as.character(l),".s",as.character(1:n))
    colnames(data.rel[[l]]$Y) <- c.name
    rownames(data.rel[[l]]$Y) <- sample.nm
  }

  ## ---- Summarize data signal
  signal.names <- colnames(data.rel[[1]]$Y)[signal.idx]

  ## ---- Batch correction by MMUPHin
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
    data.rel.mmuphin[[l]] <- list(Y = Y.pool[meta$Study == as.character(l),], X = data.rel[[l]]$X)
  }

  # ============================================================================================== #
  # Analysis on original data ----
  rm(list = setdiff(ls(), c("data.rel", "data.rel.mmuphin", "data.rel.conqur", "signal.names", "s", "data.loc", "Ka")))
  source("./utility/CRC_utility.R")
  source("./utility/roc.R")
  AUPRC = S = method = batch_type = data_type <- NULL
  L <- length(data.rel)
  K <- ncol(data.rel[[1]]$Y)
  study = Y.pool = X.pool = Y.rarify <- NULL
  for(l in 1:L){
    study <- c(study, rep(l, length(data.rel[[l]]$X)) )
    Y.pool <- rbind(Y.pool, data.rel[[l]]$Y)
    Y.rarify <- rbind(Y.rarify, .rarefy(otu.tab = data.rel[[l]]$Y, ss = 2023))
    X.pool <- c(X.pool, data.rel[[l]]$X)
  }

  ## BW (rarefied data) ----
  P.pool = Y.rarify / rowSums(Y.rarify)
  pval = NULL
  for(k in 1:K){
    BW.data = data.frame(y=P.pool[,k], x = as.factor(X.pool), block=as.factor(study) )
    pval = c(pval, pvalue(wilcox_test(y ~ x | block, data=BW.data)) )
  }

  BW.prop.model <- data.frame(p.val = pval, q.val = p.adjust(pval, method="BH"), row.names = colnames(P.pool))

  ## ---- Calculate AUPRC
  rank.fdr <- rank(BW.prop.model$p.val)
  true.active <- rep(0, nrow(BW.prop.model))
  names(true.active) <- rownames(BW.prop.model)
  true.active[signal.names] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "Original")
  data_type <- c(data_type, "Independent")
  method <- c(method, "BW")

  ## ALDEx2 ----
  ## ---- function use for ALDEx2
  source("./utility/aldex2.R")

  ## ---- Combine data
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(study))
  colnames(feature.table) <- rownames(meta.data)

  ## aldex2 will add 0.5 in it's own function, we don't need to add any pseudocount here.
  ## ---- ALDEx2 model without gamma
  Aldex2.model.default <- run_aldex2(otu.tab = feature.table, meta = meta.data, formula = "labels + study")

  ## ---- Calculate AUPRC
  true.active <- rep(0, nrow(Aldex2.model.default$glmfit))
  names(true.active) <- rownames(Aldex2.model.default$glmfit)
  true.active[signal.names] <- 1
  rank.fdr <- rank(Aldex2.model.default$glmfit$`labels1:pval`)
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "Original")
  data_type <- c(data_type, "Independent")
  method <- c(method, "ALDEx2")

  ## ANCOM-BC2 ----
  ## ---- function use for ANCOM-BC2
  source("./utility/ancombc.R")

  ## ---- ANCOM-BC2
  ANCOMBC2.model <- ancombc.fun(feature.table = feature.table,
                                meta = meta.data,
                                formula = "labels + study",
                                adjust.method = "fdr",
                                group = NULL,
                                subject = NULL,
                                method = "ancombc2")

  ## ---- Calculate AUPRC
  rank.fdr <- rank(ANCOMBC2.model$res$p_labels1)
  true.active <- rep(0, length(ANCOMBC2.model$res$taxon))
  names(true.active) <- ANCOMBC2.model$res$taxon
  true.active[intersect(signal.names, ANCOMBC2.model$res$taxon)] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "Original")
  data_type <- c(data_type, "Independent")
  method <- c(method, "ANCOM-BC2")

  ## MMUPHin ----
  ## ---- function use for MMUPHin
  source("./utility/mmuphin.R")

  ## ---- MMUPHin model
  MMUPHin.model <- fit_metaAnalysis(feature.abd = feature.table,
                                    data = meta.data,
                                    test_variable = "labels",
                                    contrasts = list("1" = "1", "0" = "0"),
                                    batch_variable = "study",
                                    covariates = NULL,
                                    covariates.random = NULL,
                                    moderator_variables = NULL)

  ## ---- Calculate AUPRC
  rank.fdr <- rank(MMUPHin.model$meta_fits$pval)
  true.active <- rep(0, nrow(MMUPHin.model$meta_fits))
  names(true.active) <- MMUPHin.model$meta_fits$feature
  true.active[signal.names] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "Original")
  data_type <- c(data_type, "Independent")
  method <- c(method, "MMUPHin")

  ## CLR-LASSO ----
  Z.pool <- log(Y.pool + 0.5)

  clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

  clrlasso <- cv.glmnet(x = clrx_Z, y = X.pool, family = 'binomial', nfolds = 5)

  ## ---- Calculate AUPRC
  true.active <- rep(0, nrow(clrlasso$glmnet.fit$beta))
  names(true.active) <- rownames(clrlasso$glmnet.fit$beta)
  true.active[signal.names] <- 1
  Paths <- list()
  for(ll in 2:ncol(clrlasso$glmnet.fit$beta)){
    Paths[[ll-1]] <- clrlasso$glmnet.fit$beta[,ll] != 0
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "Original")
  data_type <- c(data_type, "Independent")
  method <- c(method, "CLR-LASSO")

  ## Melody ----
  rel.abd <- list()
  covariate.interest <- list()
  ref <- c()
  for(l in 1:length(data.rel)){
    ref <- c(ref, "Coprococcus catus [ref_mOTU_v2_4874]")
    rel.abd[[paste0("S",l)]] <- data.rel[[l]]$Y
    covariate.interest[[paste0("S",l)]] <- data.frame(disease = data.rel[[l]]$X)
  }
  names(ref) <- names(covariate.interest)

  ## ---- Generate null model
  null.obj <- melody.null.model(rel.abd = rel.abd, ref = ref)

  ## ---- Generate summary statistics
  summary.stat.study <- melody.get.summary(null.obj = null.obj, covariate.interest = covariate.interest)

  ## ---- Meta analysis
  Melody.model <- melody.meta.summary(summary.stats = summary.stat.study,
                                      tune.path = "sequence",
                                      tune.size.sequence = 1:min((Ka + 150), 400),
                                      output.best.one = FALSE)

  ## ---- Calculate AUPRC
  true.active <- rep(0, nrow(Melody.model$disease$coef))
  names(true.active) <- rownames(Melody.model$disease$coef)
  true.active[intersect(signal.names, rownames(Melody.model$disease$coef))] <- 1
  Paths <- list()
  for(ll in 1:ncol(Melody.model$disease$coef)){
    Paths[[ll]] <- Melody.model$disease$coef[,ll]!=0
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "Original")
  data_type <- c(data_type, "Independent")
  method <- c(method, "Melody")

  # ============================================================================================== #
  # Analysis on MMUPHin batch corrected data ----
  rm(list = setdiff(ls(), c("AUPRC", "S", "data_type", "batch_type", "method", "Ka",
                            "s", "data.rel.mmuphin", "data.rel.conqur", "signal.names", "data.loc")))

  source("./utility/CRC_utility.R")
  source("./utility/roc.R")
  L <- length(data.rel.mmuphin)
  K <- ncol(data.rel.mmuphin[[1]]$Y)
  study = Y.pool = X.pool = Y.rarify <- NULL
  for(l in 1:L){
    study <- c(study, rep(l, length(data.rel.mmuphin[[l]]$X)) )
    Y.pool <- rbind(Y.pool, data.rel.mmuphin[[l]]$Y)
    Y.rarify <- rbind(Y.rarify, .rarefy(otu.tab = data.rel.mmuphin[[l]]$Y, ss = 2023))
    X.pool <- c(X.pool, data.rel.mmuphin[[l]]$X)
  }

  ## BW (rarefied data) ----
  P.pool = Y.rarify / rowSums(Y.rarify)
  pval = NULL
  for(k in 1:K){
    BW.data = data.frame(y=P.pool[,k], x = as.factor(X.pool), block=as.factor(study) )
    pval = c(pval, pvalue(wilcox_test(y ~ x | block, data=BW.data)) )
  }

  ## Multiple testing correction
  BW.prop.model <- data.frame(p.val = pval, q.val = p.adjust(pval, method="BH"), row.names = colnames(P.pool))

  ## ---- Calculate AUPRC
  rank.fdr <- rank(BW.prop.model$p.val)
  true.active <- rep(0, nrow(BW.prop.model))
  names(true.active) <- rownames(BW.prop.model)
  true.active[signal.names] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "MMUPHin")
  data_type <- c(data_type, "Independent")
  method <- c(method, "BW")

  ## ALDEx2 ----
  ## ---- function use for ALDEx2
  source("./utility/aldex2.R")

  ## ---- Combine data
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(study))
  colnames(feature.table) <- rownames(meta.data)

  ## aldex2 will add 0.5 in it's own function, we don't need to add any pseudocount here.
  ## ---- ALDEx2 model without gamma
  Aldex2.model.default <- run_aldex2(otu.tab = feature.table, meta = meta.data, formula = "labels + study")

  ## ---- Calculate AUPRC
  true.active <- rep(0, nrow(Aldex2.model.default$glmfit))
  names(true.active) <- rownames(Aldex2.model.default$glmfit)
  true.active[signal.names] <- 1
  rank.fdr <- rank(Aldex2.model.default$glmfit$`labels1:pval`)
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "MMUPHin")
  data_type <- c(data_type, "Independent")
  method <- c(method, "ALDEx2")

  ## ANCOM-BC2 ----
  ## ---- function use for ANCOM-BC2
  source("./utility/ancombc.R")

  ## ---- ANCOM-BC2
  ANCOMBC2.model <- ancombc.fun(feature.table = feature.table,
                                meta = meta.data,
                                formula = "labels + study",
                                adjust.method = "fdr",
                                group = NULL,
                                subject = NULL,
                                method = "ancombc2")

  ## ---- Calculate AUPRC
  rank.fdr <- rank(ANCOMBC2.model$res$p_labels1)
  true.active <- rep(0, length(ANCOMBC2.model$res$taxon))
  names(true.active) <- ANCOMBC2.model$res$taxon
  true.active[intersect(signal.names, ANCOMBC2.model$res$taxon)] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)

  S <- c(S, s)
  batch_type <- c(batch_type, "MMUPHin")
  data_type <- c(data_type, "Independent")
  method <- c(method, "ANCOM-BC2")

  ## MMUPHin ----
  ## ---- function use for ANCOM
  source("./utility/mmuphin.R")

  ## ---- MMUPHin default model
  MMUPHin.model <- fit_metaAnalysis(feature.abd = feature.table,
                                    data = meta.data,
                                    test_variable = "labels",
                                    contrasts = list("1" = "1", "0" = "0"),
                                    batch_variable = "study",
                                    covariates = NULL,
                                    covariates.random = NULL,
                                    moderator_variables = NULL)

  ## ---- Calculate AUPRC
  rank.fdr <- rank(MMUPHin.model$meta_fits$pval)
  true.active <- rep(0, nrow(MMUPHin.model$meta_fits))
  names(true.active) <- MMUPHin.model$meta_fits$feature
  true.active[signal.names] <- 1
  Paths <- list()
  for(ll in 1:length(rank.fdr)){
    Paths[[ll]] <- rank.fdr <= ll
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "MMUPHin")
  data_type <- c(data_type, "Independent")
  method <- c(method, "MMUPHin")

  ## CLR-LASSO ----
  Z.pool <- log(Y.pool + 0.5)

  clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

  clrlasso <- cv.glmnet(x = clrx_Z, y = X.pool, family = 'binomial', nfolds = 5)

  ## ---- Calculate AUPRC
  true.active <- rep(0, nrow(clrlasso$glmnet.fit$beta))
  names(true.active) <- rownames(clrlasso$glmnet.fit$beta)
  true.active[signal.names] <- 1
  Paths <- list()
  for(ll in 2:ncol(clrlasso$glmnet.fit$beta)){
    Paths[[ll-1]] <- clrlasso$glmnet.fit$beta[,ll] != 0
  }
  hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = FALSE)
  AUPRC <- c(AUPRC, hg.pr$AUPRC)
  S <- c(S, s)
  batch_type <- c(batch_type, "MMUPHin")
  data_type <- c(data_type, "Independent")
  method <- c(method, "CLR-LASSO")

  # ============================================================================================== #
  # Save AUPRC ----
  PRC <- data.frame(AUPRC = AUPRC, S = S, method = method, batch_type = batch_type, data_type = data_type)

  save(PRC, file = data.loc)
