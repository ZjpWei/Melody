# =============================================== #
#  CRC 2. Meta-analyses on original data (K401)   #
# =============================================== #

  # Packages ----
  library('tidyverse')
  library('miMeta')
  library('glmnet')
  library('ALDEx2')
  library("ANCOMBC")
  library('MMUPHin')
  library('coin')
  
  # General ----
  rm(list = ls())
  set.seed(2023)
  
  load("./CRC/Processed_data/data.org.K401.Rdata")
  L <- length(data.rel)
  K <- ncol(data.rel[[1]]$Y)
  study <- NULL
  Y.pool <- NULL
  X.pool <- NULL
  for(l in 1:L){
    study <- c(study, rep(l, length(data.rel[[l]]$X)) )
    Y.pool <- rbind(Y.pool, data.rel[[l]]$Y)
    X.pool <- c(X.pool, data.rel[[l]]$X)
  }
  
  # CLR-LASSO ----
  Z.pool <- log(Y.pool + 0.5)
  clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))
  clrlasso <- cv.glmnet(x = clrx_Z, y = X.pool, nfolds = 5, family = 'binomial')
  
  ## save model
  save(clrlasso, file =  "./CRC/Models_org/CLR.lasso.model.Rdata")
  
  # ALDEx2 ----
  ## functions for ALDEx2
  source("./utility/aldex2.R")

  ## combine data
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table <- data.frame(t(feature.table))
  meta.data <- data.frame(labels = as.character(outcome), study = factor(study))
  colnames(feature.table) <- rownames(meta.data)

  ## ALDEx2 model
  Aldex2.model <- run_aldex2(otu.tab = feature.table, meta = meta.data, formula = "labels + study")

  ## save model
  save(Aldex2.model, file = "./CRC/Models_org/Aldex2.model.Rdata")

  # ANCOM-BC2 ----
  ## functions for ANCOM-BC2
  source("./utility/ancombc.R")
  
  ## ANCOM-BC2 model
  ANCOMBC2.model <- ancombc.fun(feature.table = feature.table,
                                meta = meta.data,
                                formula = "labels + study",
                                adjust.method = "fdr",
                                group = "labels",
                                method = "ancombc2")
  
  ## save model
  save(ANCOMBC2.model, file = "./CRC/Models_org/ANCOMBC2.model.Rdata")
  
  # MMUPHin ----
  ## function for MMUPHin
  source("./utility/mmuphin.R")
  
  ## MMUPHin model
  MMUPHin.model <- fit_metaAnalysis(feature.abd = feature.table,
                                    data = meta.data,
                                    test_variable = "labels",
                                    contrasts = list("1" = "1", "0" = "0"),
                                    batch_variable = "study",
                                    covariates = NULL,
                                    covariates.random = NULL,
                                    moderator_variables = NULL)
  
  ## save model
  save(MMUPHin.model, file = "./CRC/Models_org/MMUPHin.model.Rdata")
  
  # BW (rarefied data) ----
  source("./utility/CRC_utility.R")
  Y.rarify <- NULL
  for(l in 1:L){
    Y.rarify <- rbind(Y.rarify, .rarefy(otu.tab = data.rel[[l]]$Y, ss = 2023))
  }

  ## Proportion normalization
  P.pool <- Y.rarify / rowSums(Y.rarify)
  pval <- NULL
  statis <- NULL
  for(k in 1:K){
    BW.data <- data.frame(y = P.pool[,k], x = factor(X.pool, levels = c("1","0")), block=as.factor(study) )
    pval <- c(pval, pvalue(wilcox_test(y ~ x | block, data = BW.data, )) )
    statis <- c(statis, statistic(wilcox_test(y ~ x | block, data = BW.data)))
  }
  
  ## multiple testing correction
  BW.prop.model <- data.frame(p.val = pval, q.val = p.adjust(pval, method="BH"),
                              statis = statis, row.names = colnames(P.pool))

  ## save model
  save(BW.prop.model, file = "./CRC/Models_org/BW.prop.model.Rdata")

  # Melody ----
  ref <- rep("Coprococcus catus [ref_mOTU_v2_4874]", L)
  rel.abd <- list()
  covariate.interest <- list()
  for(l in 1:length(data.rel)){
    rel.abd[[paste0("S",l)]] <- data.rel[[l]]$Y
    covariate.interest[[paste0("S",l)]] <- data.frame("disease" = data.rel[[l]]$X)
  }
  names(ref) <- names(rel.abd)
  
  null.obj <- melody.null.model(rel.abd = rel.abd, ref = ref)
  
  summary.stat.study <- melody.get.summary(null.obj = null.obj, covariate.interest = covariate.interest)
  
  ## meta-analysis
  Melody.model <- melody.meta.summary(summary.stats = summary.stat.study, verbose = TRUE)
  
  save(Melody.model, file = "./CRC/Models_org/Melody.model.Rdata")
  
  