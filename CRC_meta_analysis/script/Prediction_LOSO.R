# =============================================== #
#             Top X taxa & AUC lines              #
# =============================================== #

  # Packages ----
  library("tidyverse")
  library("matrixStats")
  library("paradox")
  library("lgr")
  library("progress")
  library("pROC")
  library("mlr3")
  library("mlr3learners")
  library("mlr3tuning")
  library("ranger")

  # General ----
  cat('Starting prediction script\n')
  start.time <- proc.time()[1]

  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("The analysis tag needs to be provided! Exiting...\n")
  }
  s <- as.numeric(args[1])
  ss <- as.numeric(args[2])
  methodology <- args[3]
  data.type <- args[4]

  # Parameters ----
  data.loc <- "./CRC_meta_analysis/Prediction/LOSO/"
  seed <- 2023
  train.method <- 'randomForest'
  num.folds <- 5
  num.resample <- 5
  num.trees <- c(50:100) * 10
  source("./utility/mlr3_utility.R")

  # Load data
  load(paste0(data.loc, "/prepare_data/data.rel.all.Rdata"))
  if(data.type == "Original"){
    if(!dir.exists(paste0(data.loc, "AUROC_original"))){
      dir.create(paste0(data.loc, "AUROC_original"))
    }
    load(paste0(data.loc, "index_original/index_",methodology, ".", as.character(s), ".", as.character(ss), ".Rdata"))
  }else if(data.type == "Batch-corrected"){
    if(!dir.exists(paste0(data.loc, "AUROC_batch_corrected"))){
      dir.create(paste0(data.loc, "AUROC_batch_corrected"))
    }
    load(paste0(data.loc, "index_batch_corrected/index_",methodology, ".", as.character(s), ".", as.character(ss), ".Rdata"))
  }

  # leave testing study out
  data.rel.test <- data.rel[s]
  data.rel <- data.rel[-s]
  L <- length(data.rel)
  K <- ncol(data.rel[[1]]$Y)

  # split analysis data and training data
  data.rel.analysis <- list()
  data.rel.train <- list()
  set.seed(ss + seed)
  for(l in 1:L){
    # sample case's and control's id
    case.id <- sample(which(data.rel[[l]]$X == 1))
    control.id <- sample(which(data.rel[[l]]$X == 0))
    case.rg <- round(c(0.75,1) * length(case.id))
    control.rg <- round(c(0.75,1) * length(control.id))
    
    # generate analysis data and trining data
    data.rel.analysis[[l]] <- list(Y = data.rel[[l]]$Y[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]])),],
                                   X = data.rel[[l]]$X[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]]))])

    data.rel.train[[l]] <- list(Y = data.rel[[l]]$Y[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]])),],
                                X = data.rel[[l]]$X[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]]))])
  }

  # pool training and testing data
  count.1 <- NULL
  Sample_ID <- NULL
  Study <- NULL
  Group <- NULL
  for(l in 1:length(data.rel.train)){
    count.1 <- rbind(count.1, data.rel.train[[l]]$Y)
    Sample_ID <- c(Sample_ID, rownames(data.rel.train[[l]]$Y))
    Study <- c(Study,rep(as.character(l), length(data.rel.train[[l]]$X)))
    Group <- c(Group, as.character(data.rel.train[[l]]$X))
  }
  meta.1 <- data.frame(Sample_ID = Sample_ID, Study = Study, Group = Group)
  rownames(meta.1) <- meta.1$Sample_ID

  count.2 <- NULL
  Sample_ID <- NULL
  Study <- NULL
  Group <- NULL
  for(l in 1:length(data.rel.test)){
    count.2 <- rbind(count.2, data.rel.test[[l]]$Y)
    Sample_ID <- c(Sample_ID, rownames(data.rel.test[[l]]$Y))
    Study <- c(Study,rep(as.character(l), length(data.rel.test[[l]]$X)))
    Group <- c(Group, as.character(data.rel.test[[l]]$X))
  }
  meta.2 <- data.frame(Sample_ID = Sample_ID, Study = Study, Group = Group)
  rownames(meta.2) <- meta.2$Sample_ID

  # Prediction ----
  auroc.all <- tibble(s=character(0), taxa.num=character(0),
                      method=character(0), AUC=double(0))
  block <- 'Sample_ID'
  method.names <- names(index)

  for(jj in 1:length(index)){
    j <- index[[jj]]
    # mtrys
    mtrys <- c(round(sqrt(length(j)) / 2) , round(sqrt(length(j))), round(sqrt(length(j)) * 2))

    # proportion transformation
    tmp.count.1 <- (count.1/rowSums(count.1))[,j]
    tmp.count.2 <- (count.2/rowSums(count.2))[,j]

    # split data
    set.seed(2023)
    data_splits <- create.data.split(meta = meta.1,
                                     num.folds = num.folds,
                                     num.resample = num.resample,
                                     inseparable = block)
    # train LASSO model
    model.sets <- train.model(count = tmp.count.1,
                              meta = meta.1,
                              splits = data_splits,
                              param.set=list('num.trees'= num.trees,
                                             'mtry'= mtrys)
    )

    pred_mat <- make.predictions(models = model.sets$models,
                                 count.1 = tmp.count.1,
                                 meta.1 = meta.1,
                                 count.2 = tmp.count.2)

    eval_data <- evaluate.predictions(pred_mat = pred_mat, meta = meta.2)

    auroc.all <- auroc.all %>%
      add_row(s=as.character(s),
              taxa.num=method.names[jj], method=methodology,
              AUC=eval_data$auroc %>% as.double())
  }

  if(data.type == "Original"){
    save(auroc.all, file = paste0(data.loc, "AUROC_original/auroc.",
                                  as.character(s), ".",
                                  methodology,".Rdata"))
  }else if(data.type == "Batch-corrected"){
    save(auroc.all, file = paste0(data.loc, "AUROC_batch_corrected/auroc.",
                                  as.character(s), ".",
                                  norm.method, ".",
                                  methodology,".Rdata"))
  }
