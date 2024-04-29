# =============================================== #
#             Simulation 2: Prediction            #
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
  seed <- 2024
  s <- 1
  ss <- 1
  
  # Prepare original data ----
  load("./CRC/Processed_data/data.org.K401.Rdata")
  data.rel.test <- data.rel[s]
  data.rel.loso <- data.rel[-s]
  L <- length(data.rel.loso)
  K <- ncol(data.rel.loso[[1]]$Y)
  
  data.rel.analysis <- list()
  data.rel.train <- list()
  split.info <- list()
  set.seed(ss + seed)
  for(l in 1:L){
    case.id <- sample(which(data.rel.loso[[l]]$X == 1))
    control.id <- sample(which(data.rel.loso[[l]]$X == 0))
    split.info[[l]] <- list(case.id = case.id, control.id = control.id)
    
    case.rg <- round(c(0.75,1) * length(case.id))
    control.rg <- round(c(0.75,1) * length(control.id))
    data.rel.analysis[[l]] <- list(Y = data.rel.loso[[l]]$Y[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]])),],
                                   X = data.rel.loso[[l]]$X[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]]))])
    
    data.rel.train[[l]] <- list(Y = data.rel.loso[[l]]$Y[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]])),],
                                X = data.rel.loso[[l]]$X[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]]))])
  }
  
  # Original data analysis ----
  Org.indices <- list()
  
  source("./utility/rarify.R")
  study <- NULL
  Y.pool <- NULL
  X.pool <- NULL
  Y.rarify <- NULL
  rel.abd <- list()
  covariate.interest <- list()
  for(l in 1:L){
    study <- c(study, rep(l, length(data.rel.analysis[[l]]$X)) )
    Y.pool <- rbind(Y.pool, data.rel.analysis[[l]]$Y)
    Y.rarify <- rbind(Y.rarify, .rarefy(otu.tab = data.rel.analysis[[l]]$Y, ss = 2023))
    X.pool <- c(X.pool, data.rel.analysis[[l]]$X)
    rel.abd[[paste0("S",l)]] <- data.rel.analysis[[l]]$Y
    covariate.interest[[paste0("S",l)]] <- data.frame("disease" = data.rel.analysis[[l]]$X)
  }

  ## ALDEx2
  ## function use for ALDEx2
  source("./utility/aldex2.R")
  
  ## combine data
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(study))
  colnames(feature.table) <- rownames(meta.data)
  
  ## aldex2 will add 0.5 in it's own function, we don't need to add any pseudocount here.
  Aldex2.model <- run_aldex2(otu.tab = feature.table, meta = meta.data, formula = "labels + study")
  
  index <- list()
  for(ll in c(3:16) * 5){
    index[[as.character(ll)]] <- rownames(Aldex2.model$glmfit)[rank(Aldex2.model$glmfit$`labels1:pval`) <= ll]
  }
  Org.indices[["ALDEx2"]] <- index
  
  ## ANCOM-BC2
  ## function use for ANCOMBC
  source("./utility/ancombc.R")
  
  ## ANCOM-BC2
  ANCOMBC2.model <- ancombc.fun(feature.table = feature.table,
                                meta = meta.data,
                                formula = "labels + study",
                                adjust.method = "fdr",
                                group = NULL,
                                subject = NULL,
                                method = "ancombc2")
  
  index <- list()
  for(ll in c(3:16) * 5){
    index[[as.character(ll)]] <- ANCOMBC2.model$res$taxon[rank(ANCOMBC2.model$res$p_labels1) <= ll]
  }
  Org.indices[["ANCOM-BC2"]] <- index

  ## MMUPHin
  ## function use for MMUPHin
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
  
  index <- list()
  for(ll in c(3:16) * 5){
    index[[as.character(ll)]] <- MMUPHin.model$meta_fits$feature[rank(MMUPHin.model$meta_fits$pval) <= ll]
  }
  Org.indices[["MMUPHin"]] <- index
 
  ## CLR-LASSO
  Z.pool <- log(Y.pool + 0.5)
  
  clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))
  
  clrlasso <- cv.glmnet(x = clrx_Z, y = X.pool, family = 'binomial', nfolds = 5)
  
  lentmp <- colSums(clrlasso$glmnet.fit$beta!=0)
  index <- list()
  for(ll in c(3:16) * 5){
    tmp.id <- rank(-abs(clrlasso$glmnet.fit$beta[,min(which(lentmp >= ll))])) <= ll
    index[[as.character(ll)]] <- rownames(clrlasso$glmnet.fit$beta)[tmp.id]
  }
  Org.indices[["CLR-LASSO"]] <- index
  
  ## BW (rarefied data)
  ## Log-proportion
  P.pool <- Y.rarify / rowSums(Y.rarify)
  pval <- NULL
  for(k in 1:K){
    BW.data = data.frame(y=P.pool[,k], x = as.factor(X.pool), block=as.factor(study) )
    pval = c(pval, pvalue(wilcox_test(y ~ x | block, data=BW.data)) )
  }
  
  ## Multiple testing correction
  BW.prop.model <- data.frame(p.val = pval, q.val = p.adjust(pval, method="BH"),
                              row.names = colnames(P.pool))
  
  index <- list()
  for(ll in c(3:16) * 5){
    index[[as.character(ll)]] <- rownames(BW.prop.model)[rank(BW.prop.model$p.val) <= ll]
  }
  Org.indices[["BW"]] <- index
  
  ## Melody
  ## null model
  null.obj <- melody.null.model(rel.abd = rel.abd)
  
  ## summary statistics
  summary.stat.study <- melody.get.summary(null.obj = null.obj, covariate.interest = covariate.interest)
  
  ## meta-analysis
  Melody.model <- melody.meta.summary(summary.stats = summary.stat.study, 
                                      tune.path = "sequence",
                                      tune.size.sequence = 1:80,
                                      output.best.one = FALSE)
  
  index <- list()
  for(ll in c(3:16) * 5){
    index[[as.character(ll)]] <- rownames(Melody.model$disease$coef)[Melody.model$disease$coef[,ll]!=0]
  }
  Org.indices[["Melody"]] <- index
  
  # Make predictions
  source("./utility/Make_prediction.R")
  AUROC.org <- tibble(s=character(0), taxa.num=character(0), method=character(0), AUC=double(0))
  
  for(nm.id in names(Org.indices)){
    AUROC.org <- AUROC.org %>% add_row(
      Make_prediction(data.rel.train = data.rel.train, data.rel.test = data.rel.test, 
                      index = Org.indices[[nm.id]], method = nm.id))
  }
  rm(list = setdiff(ls(), c("s", "ss", "split.info", "AUROC.org")))
  
  # Prepare batch-corrected data ----
  load("./CRC/Processed_data/data.bat.K401.Rdata")
  data.rel.test <- data.rel.batch[s]
  data.rel.loso <- data.rel.batch[-s]
  L <- length(data.rel.loso)
  K <- ncol(data.rel.loso[[1]]$Y)
  
  data.rel.analysis <- list()
  data.rel.train <- list()
  for(l in 1:L){
    case.id <- split.info[[l]]$case.id
    control.id <- split.info[[l]]$control.id
    
    case.rg <- round(c(0.75,1) * length(case.id))
    control.rg <- round(c(0.75,1) * length(control.id))
    data.rel.analysis[[l]] <- list(Y = data.rel.loso[[l]]$Y[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]])),],
                                   X = data.rel.loso[[l]]$X[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]]))])
    
    data.rel.train[[l]] <- list(Y = data.rel.loso[[l]]$Y[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]])),],
                                X = data.rel.loso[[l]]$X[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]]))])
  }
  
  # Original data analysis ----
  Bat.indices <- list()
  
  source("./utility/rarify.R")
  study <- NULL
  Y.pool <- NULL
  X.pool <- NULL
  Y.rarify <- NULL
  for(l in 1:L){
    study <- c(study, rep(l, length(data.rel.analysis[[l]]$X)) )
    Y.pool <- rbind(Y.pool, data.rel.analysis[[l]]$Y)
    Y.rarify <- rbind(Y.rarify, .rarefy(otu.tab = data.rel.analysis[[l]]$Y, ss = 2023))
    X.pool <- c(X.pool, data.rel.analysis[[l]]$X)
  }
  
  ## ALDEx2
  ## function use for ALDEx2
  source("./utility/aldex2.R")
  
  ## combine data
  outcome <- X.pool
  feature.table <- Y.pool
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(study))
  colnames(feature.table) <- rownames(meta.data)
  
  ## aldex2 will add 0.5 in it's own function, we don't need to add any pseudocount here.
  Aldex2.model <- run_aldex2(otu.tab = feature.table, meta = meta.data, formula = "labels + study")
  
  index <- list()
  for(ll in c(3:16) * 5){
    index[[as.character(ll)]] <- rownames(Aldex2.model$glmfit)[rank(Aldex2.model$glmfit$`labels1:pval`) <= ll]
  }
  Bat.indices[["ALDEx2"]] <- index
  
  ## ANCOM-BC2
  ## function use for ANCOMBC
  source("./utility/ancombc.R")
  
  ## ANCOM-BC2
  ANCOMBC2.model <- ancombc.fun(feature.table = feature.table,
                                meta = meta.data,
                                formula = "labels + study",
                                adjust.method = "fdr",
                                group = NULL,
                                subject = NULL,
                                method = "ancombc2")
  
  index <- list()
  for(ll in c(3:16) * 5){
    index[[as.character(ll)]] <- ANCOMBC2.model$res$taxon[rank(ANCOMBC2.model$res$p_labels1) <= ll]
  }
  Bat.indices[["ANCOM-BC2"]] <- index
  
  ## MMUPHin
  ## function use for MMUPHin
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
  
  index <- list()
  for(ll in c(3:16) * 5){
    index[[as.character(ll)]] <- MMUPHin.model$meta_fits$feature[rank(MMUPHin.model$meta_fits$pval) <= ll]
  }
  Bat.indices[["MMUPHin"]] <- index
  
  ## CLR-LASSO
  Z.pool <- log(Y.pool + 0.5)
  
  clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))
  
  clrlasso <- cv.glmnet(x = clrx_Z, y = X.pool, family = 'binomial', nfolds = 5)
  
  lentmp <- colSums(clrlasso$glmnet.fit$beta!=0)
  index <- list()
  for(ll in c(3:16) * 5){
    tmp.id <- rank(-abs(clrlasso$glmnet.fit$beta[,min(which(lentmp >= ll))])) <= ll
    index[[as.character(ll)]] <- rownames(clrlasso$glmnet.fit$beta)[tmp.id]
  }
  Bat.indices[["CLR-LASSO"]] <- index
  
  ## BW (rarefied data)
  ## Log-proportion
  P.pool <- Y.rarify / rowSums(Y.rarify)
  pval <- NULL
  for(k in 1:K){
    BW.data = data.frame(y=P.pool[,k], x = as.factor(X.pool), block=as.factor(study) )
    pval = c(pval, pvalue(wilcox_test(y ~ x | block, data=BW.data)) )
  }
  
  ## Multiple testing correction
  BW.prop.model <- data.frame(p.val = pval, q.val = p.adjust(pval, method="BH"),
                              row.names = colnames(P.pool))
  
  index <- list()
  for(ll in c(3:16) * 5){
    index[[as.character(ll)]] <- rownames(BW.prop.model)[rank(BW.prop.model$p.val) <= ll]
  }
  Bat.indices[["BW"]] <- index
  
  # Make predictions
  source("./utility/Make_prediction.R")
  AUROC.bat <- tibble(s=character(0), taxa.num=character(0), method=character(0), AUC=double(0))
  
  for(nm.id in names(Bat.indices)){
    AUROC.bat <- AUROC.bat %>% add_row(
      Make_prediction(data.rel.train = data.rel.train, data.rel.test = data.rel.test, 
                      index = Bat.indices[[nm.id]], method = nm.id))
  }
  rm(list = setdiff(ls(), c("s", "ss", "split.info", "AUROC.org",  "AUROC.bat")))
  
  # save data ----
  AUROC.org$type <- "Original"
  AUROC.bat$type <- "Study-effect corrected"
  AUROC.all <- AUROC.org %>% add_row(AUROC.bat)
  
  save(AUROC.all, file = paste0("./AUROC/auroc.", as.character(s), ".",as.character(ss), ".Rdata"))
