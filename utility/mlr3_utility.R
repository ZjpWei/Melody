create.data.split <- function(meta,
                              num.folds = 5, 
                              num.resample = 5,
                              type = "BINARY",
                              stratify = TRUE, 
                              inseparable = NULL, 
                              verbose = 1) {
  ####
  tmp.lb <- as.numeric(meta$Group) * 2 - 1
  names(tmp.lb) <- rownames(meta)
  tmp.info <- c(-1,1)
  names(tmp.info) <- as.character(c(0,1))
  ####
  label <- list(label = tmp.lb,
                info = tmp.info,
                type = type)
  
  group.numbers <- vapply(label$info,
                          FUN = function(x){
                            sum(label$label == x)},
                          FUN.VALUE = integer(1))
  if (any(group.numbers <= 5)){
    msg <- paste0("Data set has only:\n",
                  paste0(names(group.numbers)[1], "\t", group.numbers[1]),
                  "\n",
                  paste0(names(group.numbers)[2], "\t", group.numbers[2]),
                  "\nThis is not enough for SIAMCAT to proceed!")
    stop(msg)
  }
  
  labelNum <- as.numeric(label$label)
  names(labelNum) <- names(label$label)
  exm.ids <- names(labelNum)
  
  if (is.null(inseparable) || inseparable == "" ||
      toupper(inseparable) == "NULL" ||
      toupper(inseparable) == "NONE" ||
      toupper(inseparable) == "UNKNOWN") {
    inseparable <- NULL
  }
  
  # parse label description
  classes <- sort(label$info)
  
  ### check arguments
  if (num.resample < 1) {
    if (verbose > 1){
      msg <- paste0("+++ Resetting num.resample = 1 (", 
                    num.resample, 
                    " is an invalid number of resampling rounds)")
      message(msg)
    }
    num.resample <- 1
  }
  if (num.folds < 2) {
    if (verbose > 1){
      msg <- paste0("+++ Resetting num.folds = 2 (", 
                    num.folds, " is an invalid number of folds)")
      message(msg)
    }
    num.folds <- 2
  }
  if (!is.null(inseparable) && stratify) {
    if (verbose > 1){
      msg <- paste0("+++ Resetting stratify to FALSE ",
                    "(Stratification is not supported when ",
                    "inseparable is given")
      message(msg)
    }
    stratify <- FALSE
  }
  if (num.folds >= length(labelNum)) {
    if (verbose > 1)
      message("+++ Performing un-stratified",
              "leave-one-out (LOO) cross-validation")
    stratify <- FALSE
    num.folds <- length(labelNum)
  }
  if (!is.null(inseparable) && is.null(meta)) {
    stop("Meta-data must be provided if the inseparable parameter is not
                NULL")
  }
  if (!is.null(inseparable)) {
    if (is.numeric(inseparable) && length(inseparable) == 1) {
      stopifnot(inseparable <= ncol(meta))
    } else if (is.character(inseparable) &&
               length(inseparable == 1)) {
      stopifnot(inseparable %in% colnames(meta))
    } else {
      stop(
        "Inseparable parameter must be either a single column index
                    or a single column name of metadata matrix"
      )
    }
  }
  
  train.list <- list(NULL)
  test.list <- list(NULL)
  
  
  for (r in seq_len(num.resample)) {
    labelNum <- sample(labelNum)
    if (label$type == 'BINARY'){
      foldid <-
        assign.fold.binary(
          label = labelNum,
          num.folds = num.folds,
          stratified = stratify,
          inseparable = inseparable,
          meta = meta[names(labelNum),],
          verbose = verbose)
    } else if (label$type == 'CONTINUOUS'){
      foldid <-
        assign.fold.regr(
          label = labelNum,
          num.folds = num.folds,
          inseparable = inseparable,
          meta = meta[names(labelNum),],
          verbose = verbose)
    }
    names(foldid) <- names(labelNum)
    stopifnot(length(labelNum) == length(foldid))
    stopifnot(length(unique(foldid)) == num.folds)
    
    train.temp <- list(NULL)
    test.temp <- list(NULL)
    
    if (verbose > 1){
      msg <- paste("+ resampling round", r)
      message(msg)
    }
    for (f in seq_len(num.folds)) {
      # make sure each fold contains examples from all classes for
      # stratify==TRUE should be tested before assignment of
      # test/training set
      if (stratify) {
        stopifnot(all(sort(unique(labelNum[foldid == f])) ==
                        classes))
      }
      # select test examples
      test.idx <- which(foldid == f)
      train.idx <- which(foldid != f)
      train.temp[f] <- list(names(foldid)[train.idx])
      test.temp[f] <- list(names(foldid)[test.idx])
      # for startify==FALSE, all classes must only be present in the
      # training set e.g. in leave-one-out CV, the test fold
      # cannot contain all classes
      if (!stratify && label$type == 'BINARY') {
        stopifnot(all(sort(unique(labelNum[foldid != f]))
                      == classes))
      }
      stopifnot(length(intersect(train.idx, test.idx)) == 0)
      if (verbose > 2){
        msg <- paste("+++ fold ", f, " contains ",
                     sum(foldid == f), " samples")
        message(msg)
      }
    }
    train.list[[r]] <- train.temp
    test.list[[r]] <- test.temp
  }
  
  data_split <- list(
    training.folds = train.list,
    test.folds = test.list,
    num.resample = num.resample,
    num.folds = num.folds
  )
  e.time <- proc.time()[3]
  return(data_split)
}


#' @keywords internal
assign.fold.binary <- function(label, num.folds, stratified,
                               inseparable = NULL, meta = NULL, verbose = 1) {
  if (verbose > 2)
    message("+++ starting assign.fold.binary")
  foldid <- rep(0, length(label))
  classes <- sort(unique(label))
  # Transform number of classes into vector of 1 to x for looping over.
  # stratify positive examples
  if (stratified) {
    # If stratify is TRUE, make sure that num.folds does not exceed the
    # maximum number of examples for the class with
    # the fewest training examples.
    if (any(as.data.frame(table(label))[, 2] < num.folds)) {
      stop(
        "+++ Number of CV folds is too large for this data set to
                    maintain stratification. Reduce num.folds or turn
                    stratification off. Exiting."
      )
    }
    for (c in seq_along(classes)) {
      idx <- which(label == classes[c])
      foldid[idx] <- sample(rep(seq_len(num.folds),
                                length.out = length(idx)))
    }
  } else {
    # If stratify is not TRUE, make sure that num.sample is not
    # bigger than number.folds
    if (length(label) < num.folds) {
      warning(
        "+++ num.samples is exceeding number of folds,",
        " setting CV to (k-1) unstratified CV"
      )
      num.folds <- length(label)
    }
    if (!is.null(inseparable)) {
      strata <- unique(meta[[inseparable]])
      sid <-
        sample(rep(seq_len(num.folds), length.out =
                     length(strata)))
      for (s in seq_along(strata)) {
        idx <- which(meta[[inseparable]] == strata[s])
        foldid[idx] <- sid[s]
      }
      stopifnot(all(!is.na(foldid)))
    } else {
      foldid <- sample(rep(seq_len(num.folds),
                           length.out = length(label)))
    }
  }
  # make sure that for each test fold the training fold (i.e. all other
  # folds together) contain examples from all classes except for
  # stratified CV
  if (!stratified) {
    for (f in seq_len(num.folds)) {
      stopifnot(all(sort(unique(label[foldid != f])) == classes))
    }
  } else {
    for (f in seq_len(num.folds)) {
      stopifnot(all(sort(unique(label[foldid == f])) == classes))
    }
  }
  
  stopifnot(length(label) == length(foldid))
  if (verbose > 2)
    message("+++ finished assign.fold.binary")
  return(foldid)
}


#' @keywords internal
assign.fold.regr <- function(label, num.folds, inseparable = NULL,
                             meta = NULL, verbose = 1) {
  if (verbose > 2)
    message("+++ starting assign.fold.regr")
  foldid <- rep(0, length(label))
  
  # If stratify is not TRUE, make sure that num.sample is not
  # bigger than number.folds
  
  if (!is.null(inseparable)) {
    strata <- unique(meta[[inseparable]])
    sid <- sample(rep(seq_len(num.folds), length.out = length(strata)))
    for (s in seq_along(strata)) {
      idx <- which(meta[[inseparable]] == strata[s])
      foldid[idx] <- sid[s]
    }
    stopifnot(all(!is.na(foldid)))
  } else {
    foldid <- sample(rep(seq_len(num.folds), length.out = length(label)))
  }
  
  stopifnot(length(label) == length(foldid))
  if (verbose > 2)
    message("+++ finished assign.fold.regr")
  return(foldid)
}


train.model <- function(count,
                        meta,
                        splits,
                        type = "BINARY",
                        method = "randomForest",
                        measure = "classif.acc", 
                        param.set = NULL,
                        grid.size=11, 
                        feature.type='normalized', 
                        verbose = 1) {
  
  # check and get features
  feat <- t(count)
  # make sure the names fit
  rownames(feat) <- make.names(rownames(feat))
  # checks
  ####
  tmp.lb <- as.numeric(meta$Group) * 2 - 1
  names(tmp.lb) <- rownames(meta)
  tmp.info <- c(-1,1)
  names(tmp.info) <- as.character(c(0,1))
  ####
  label <- list(label = tmp.lb,
                info = tmp.info,
                type = type)
  
  data.split <- splits
  # Create List to save models.
  models.list <- list()
  num.runs <- data.split$num.folds * data.split$num.resample
  
  ## Create Learner
  lrn <- create.mlr.learner(method, nrow(feat), 
                            param.set,
                            type=label$type)
  
  
  get_logger("mlr3")$set_threshold('off')
  get_logger("bbotk")$set_threshold('off')
  # loop over the folds
  bar <- 0
  if (verbose > 1){
    msg <- paste("+ training", method, "models on", num.runs,
                 "training sets")
    message(msg)
  }
  
  if (verbose == 1 || verbose == 2)
    pb <- progress_bar$new(total = num.runs)
  for (fold in seq_len(data.split$num.folds)) {
    if (verbose > 2){
      msg <- paste("+++ training on cv fold:", fold)
      message(msg)
    }
    
    for (resampling in seq_len(data.split$num.resample)) {
      if (verbose > 2){
        msg <- paste("++++ repetition:", resampling)
        message(msg)
      }
      ## Prepare data
      fold.name <- paste0("cv_fold", as.character(fold), "_rep",
                          as.character(resampling))
      fold.exm.idx <-
        match(data.split$training.folds[[resampling]][[fold]],
              names(label$label))
      
      ### subselect examples for training
      label.fac <- label$label
      if (label$type == 'BINARY'){
        label.fac <- factor(label.fac,  levels = sort(label$info))
      }
      train.label <- label.fac[fold.exm.idx]
      data <-
        as.data.frame(t(feat)[fold.exm.idx,])
      stopifnot(nrow(data) == length(train.label))
      stopifnot(all(rownames(data) == names(train.label)))
      
      data$label <- train.label
      
      if (label$type == 'BINARY' && data$label[1] != -1) {
        data <- data[c(which(data$label == -1)[1],
                       c(seq_len(nrow(data)))[-which(data$label == -1)[1]]), ]
      }
      
      ## create task
      if (label$type == 'BINARY'){
        task <- TaskClassif$new(id='classif', backend=data,
                                target='label', positive="1")
      } else if (label$type == 'CONTINUOUS') {
        task <- TaskRegr$new(id='regr', backend=data, target='label')
      }
      lrn.fold <- lrn$clone(deep=TRUE)
      
      ## Train model
      any.tuner <- unlist(lapply(lrn$param_set$values, FUN=class))
      if (any(any.tuner=='TuneToken')){
        instance <- tune(tnr("grid_search", resolution=grid.size),
                         task = task,
                         learner = lrn.fold,
                         resampling = rsmp("cv", folds = 5),
                         measures = msr(measure))
        lrn.fold$param_set$values <- instance$result_learner_param_vals
      }
      model <- lrn.fold$train(task = task)
      feat.weights <- model$importance()
      # add feature weights to the model
      bar <- bar + 1
      models.list[[fold.name]] <- list(model=model, features=feat.weights)
      if (verbose == 1 || verbose == 2)
        pb$tick()
    }
  }
  
  model_list <- list(
    models = models.list,
    model.type = method,
    feature.type = feature.type)
  
  e.time <- proc.time()[3]
  
  if (verbose > 1){
    msg <- paste("+ finished train.model in", 
                 formatC(e.time - s.time, digits = 3), "s")
    message(msg)
  }
  if (verbose == 1){
    msg <- paste("Trained", method, "models successfully.")
    message(msg)
  }
  
  return(model_list)
}

make.predictions <- function(models, 
                             count.1,
                             meta.1,
                             count.2,
                             type = "BINARY",
                             verbose = 1) {
  if (verbose > 1)
    message("+ starting make.predictions on external dataset")
  
  #models <- models(siamcat.trained)
  ####
  tmp.lb <- as.numeric(meta.1$Group) * 2 - 1
  names(tmp.lb) <- rownames(meta.1)
  tmp.info <- unique(tmp.lb)
  names(tmp.info) <- as.character((tmp.info + 1) / 2)
  ####
  label <- list(label = tmp.lb,
                info = tmp.info,
                type = type)
  
  feat.test <- count.2
  colnames(feat.test) <- make.names(colnames(feat.test))
  
  feat.ref <- count.1
  colnames(feat.ref) <- make.names(colnames(feat.ref))
  
  # data sanity checks
  stopifnot(all(colnames(feat.ref) %in% colnames(feat.test)))
  # prediction
  num.models <- length(models)
  
  pred <- matrix(NA, ncol = num.models, nrow = nrow(feat.test),
                 dimnames = list(rownames(feat.test),
                                 paste0("Model_", seq_len(num.models))))
  if (verbose == 1 || verbose == 2)
    pb <- progress_bar$new(
      total = num.models)
  for (i in seq_len(num.models)) {
    data <- as.data.frame(feat.test)
    model <- models[[i]]
    
    data <- data[, names(model$features)]
    if (verbose > 2){
      msg <- paste0("Applying ", model_type(siamcat.trained),
                    " on complete external dataset", " (", i, " of ",
                    num.models, ")...")
      message(msg)
    }
    if (label$type == 'BINARY'){
      data$label <- as.factor(c(unique(label$label),
                                sample(unique(label$label), 
                                       size=nrow(data)-2, replace=TRUE)))
      test.task <- TaskClassif$new(id='classif', backend=data,
                                   target='label')
      pdata <- model$model$predict(task=test.task)
      p <- pdata$data$prob[,2]
      names(p) <- rownames(data)
    } else if (label$type == 'CONTINUOUS'){
      data$label <- 0
      test.task <- TaskRegr$new(id='regr', backend=data, target='label')
      pdata <- model$model$predict(task=test.task)
      p <- pdata$data$response
      names(p) <- rownames(data)
    }
    
    pred[names(p), i] <- p
    
    if (verbose == 1 || verbose == 2)
      pb$tick()
  }
  pred_matrix <- pred
  return(pred_matrix)
}

evaluate.predictions <- function(pred_mat,
                                 meta, 
                                 type = "BINARY",
                                 verbose = 1) {
  
  s.time <- proc.time()[3]
  ####
  tmp.lb <- as.numeric(meta$Group) * 2 - 1
  names(tmp.lb) <- rownames(meta)
  tmp.info <- c(-1,1)
  names(tmp.info) <- as.character(c(0,1))
  ####
  label <- list(label = tmp.lb,
                info = tmp.info,
                type = type)
  
  if (label$type == 'BINARY'){
    r.object <- eval.binary(pred_mat = pred_mat, label = label, s.time, verbose)
  } else if (label$type == 'CONTINUOUS'){
    r.object <- eval.regr(pred_mat = pred_mat, label = label, s.time, verbose)
  }
  
  return(r.object)
}


#' @keywords internal
eval.regr <- function(pred_mat,
                      label,
                      s.time, 
                      verbose=0){
  
  m <- match(names(label$label), rownames(pred_mat))
  
  pred <- pred_mat[m, , drop = FALSE]
  stopifnot(all(names(label$label) == rownames(pred)))
  
  if (ncol(pred) > 1){
    if (verbose > 2)
      message("+ evaluating multiple predictions")
    
    # mean predictions
    pred.mean <- rowMeans(pred)
    ess <- sum((label$label - mean(label$label))^2)
    rss <- sum((abs(pred.mean) - label$label)^2)
    r2.mean <- 1 - rss/ess
    mae.mean <- mean(abs(pred.mean - label$label))
    mse.mean <- mean((abs(pred.mean - label$label))^2)
    
    # each repetition individually
    r2.all <- list()
    mae.all <- list()
    mse.all <- list()
    for (i in seq_len(ncol(pred))){
      rss <- sum((abs(pred[,i]) - label$label)^2)
      r2.all[[i]] <- 1 - rss/ess
      mae.all[[i]] <- mean(abs(pred[,i] - label$label))
      mse.all[[i]] <- mean((abs(pred[,i] - label$label))^2)
    }
    
    eval_data <- list(
      r2 = r2.mean, r2.all = r2.all,
      mae = mae.mean, mae.all = mae.all,
      mse = mse.mean, mse.all = mse.all)
    
  } else {
    if (verbose > 2)
      message("+ evaluating single prediction")
    rss <- sum((abs(pred[,1]) - label$label)^2)
    ess <- sum((label$label - mean(label$label))^2)
    r2 <- 1 - rss/ess
    mae <- mean(abs(pred[,1] - label$label))
    mse <- mean((abs(pred[,1] - label$label))^2)
    eval_data <- list(r2=r2, mae=mae, mse=mse)
  }
  e.time <- proc.time()[3]
  if (verbose > 1){
    msg <- paste(
      "+ finished evaluate.predictions in",
      formatC(e.time - s.time, digits = 3),
      "s")
    message(msg)
  }
  if (verbose == 1)
    message("Evaluated predictions successfully.")
  return(eval_data)
}

eval.binary <- function(pred_mat,
                        label,
                        s.time, 
                        verbose=0){
  
  summ.stat <- "mean" # TODO make this a possible parameter?
  # TODO compare header to label make sure that label and prediction are in
  # the same order
  m <- match(names(label$label), rownames(pred_mat))
  
  pred <- pred_mat[m, , drop = FALSE]
  stopifnot(all(names(label$label) == rownames(pred)))
  
  # ##########################################################################
  # ROC curve
  if (verbose > 2)
    message("+ calculating ROC")
  auroc <- 0
  if (ncol(pred) > 1) {
    roc.all <- list()
    auroc.all <- vector("numeric", ncol(pred))
    for (c in seq_len(ncol(pred))) {
      roc.all[[c]] <- roc(response = label$label, predictor = pred[, c],
                          direction = '<', levels = label$info, ci = FALSE)
      auroc.all[c] <- roc.all[[c]]$auc
    }
    l.vec <- rep(label$label, ncol(pred))
  } else {
    l.vec <- label$label
  }
  
  # average data for plotting one mean prediction curve
  
  roc.mean <- roc(response = label$label,
                  predictor = apply(pred, 1, summ.stat),
                  ci = TRUE, of = "se",
                  sp = seq(0, 1, 0.05), direction = '<', levels = label$info)
  auroc <- roc.mean$auc
  
  # ##########################################################################
  # PR curve
  prc <- list()
  ev <- list()
  auprc <- 0
  if (ncol(pred) > 1) {
    auprc.all <- vector("numeric", ncol(pred))
    prc.all <- list()
    ev.all <- list()
    for (c in seq_len(ncol(pred))) {
      ev.all[[c]] <- evaluate.classifier(pred[, c], label$label, label,
                                         verbose = verbose)
      prc.all[[c]] <- evaluate.get.pr(ev.all[[c]], verbose = verbose)
      auprc.all[c] <- evaluate.calc.aupr(ev.all[[c]], verbose = verbose)
    }
    ev <- evaluate.classifier(apply(pred, 1, summ.stat), label$label, label)
  } else {
    ev <- evaluate.classifier(as.vector(pred), label$label, label,
                              verbose = verbose)
  }
  
  prc <- evaluate.get.pr(ev, verbose = verbose)
  auprc <- c(evaluate.calc.aupr(ev, verbose = verbose))
  
  if (ncol(pred) > 1) {
    if (verbose > 2)
      message("+ evaluating multiple predictions")
    eval_data <- list(
      roc= roc.mean, roc.all = roc.all,
      auroc = auroc, auroc.all = auroc.all,
      prc = prc, prc.all = prc.all,
      auprc = auprc, auprc.all = auprc.all,
      ev = ev, ev.all = ev.all
    )
    
  } else {
    if (verbose > 2)
      message("+ evaluating single prediction")
    eval_data <- list(
      roc=roc.mean, auroc=auroc, prc=prc, auprc=auprc, ev=ev
    )
  }
  e.time <- proc.time()[3]
  if (verbose > 1){
    msg <- paste(
      "+ finished evaluate.predictions in",
      formatC(e.time - s.time, digits = 3),
      "s")
    message(msg)
  }
  if (verbose == 1)
    message("Evaluated predictions successfully.")
  return(eval_data)
}

evaluate.classifier <-
  function(predictions, test.label, label, verbose = 0) {
    if (verbose > 2)
      message("+ starting evaluate.classifier")
    stopifnot(is.null(dim(test.label)))
    stopifnot(length(unique(test.label)) == 2)
    stopifnot(all(is.finite(predictions)))
    # calculate thresholds, one between each subsequent pair of sorted
    #prediction values this is ignorant to whether
    # predictions is in matrix or vector format (see below)
    thr <- predictions
    dim(thr) <- NULL
    thr <- sort(unique(thr))
    thr <-
      rev(c(min(thr) - 1, (thr[-1] + thr[-length(thr)]) / 2,
            max(thr) + 1))
    if (is.null(dim(predictions))) {
      # assuming that a single model was applied to predict the data set
      stopifnot(length(test.label) == length(predictions))
      # actual evaluations per threshold value
      tp <- vapply(
        thr,
        FUN = function(x) {
          sum(test.label == max(label$info)
              & predictions > x)
        },
        USE.NAMES = FALSE,
        FUN.VALUE = integer(1)
      )
      fp <- vapply(
        thr,
        FUN = function(x) {
          sum(test.label == min(label$info)
              & predictions > x)
        },
        USE.NAMES = FALSE,
        FUN.VALUE = integer(1)
      )
      tn <- vapply(
        thr,
        FUN = function(x) {
          sum(test.label == min(label$info)
              & predictions < x)
        },
        USE.NAMES = FALSE,
        FUN.VALUE = integer(1)
      )
      fn <- vapply(
        thr,
        FUN = function(x) {
          sum(test.label == max(label$info)
              & predictions < x)
        },
        USE.NAMES = FALSE,
        FUN.VALUE = integer(1)
      )
    } else {
      # assuming that several models were applied to predict the same data
      # and predictions of each model occupy one
      # column
      stopifnot(length(test.label) == nrow(predictions))
      tp <- t(vapply(
        thr,
        FUN = function(x) {
          apply(
            predictions,
            2,
            FUN = function(y) {
              sum(test.label == max(label$info) & y > x)
            }
          )
        },
        USE.NAMES = FALSE,
        FUN.VALUE = integer(2)
      ))
      fp <- t(vapply(
        thr,
        FUN = function(x) {
          apply(
            predictions,
            2,
            FUN = function(y) {
              sum(test.label == min(label$info) & y > x)
            }
          )
        },
        USE.NAMES = FALSE,
        FUN.VALUE = integer(2)
      ))
      tn <- t(vapply(
        thr,
        FUN = function(x) {
          apply(
            predictions,
            2,
            FUN = function(y) {
              sum(test.label == min(label$info) & y < x)
            }
          )
        },
        USE.NAMES = FALSE,
        FUN.VALUE = integer(2)
      ))
      fn <- t(vapply(
        thr,
        FUN = function(x) {
          apply(
            predictions,
            2,
            FUN = function(y) {
              sum(test.label == max(label$info) & y < x)
            }
          )
        },
        USE.NAMES = FALSE,
        FUN.VALUE = integer(2)
      ))
    }
    if (verbose > 2)
      message("+ finished evaluate.classifier")
    return(list(
      tp = tp,
      tn = tn,
      fp = fp,
      fn = fn,
      thresholds = thr
    ))
  }

# calculates the area under a curve using a trapezoid approximation
evaluate.area.trapez <- function(x, y, verbose = 0) {
  if (verbose > 2)
    message("+ starting evaluate.area.trapez")
  if (x[1] > x[length(x)]) {
    x <- rev(x)
    y <- rev(y)
  }
  xd <- x[-1] - x[-length(x)]
  ym <- 0.5 * (y[-1] + y[-length(y)])
  if (verbose > 2)
    message("+ finished evaluate.area.trapez")
  return(xd %*% ym)
}

evaluate.get.pr <- function(eval, verbose = 0) {
  if (verbose > 2)
    message("+ starting evaluate.get.pr")
  tpr <- eval$tp / (eval$tp + eval$fn)
  ppv <- eval$tp / (eval$tp + eval$fp)
  # at thresholds where the classifier makes no positive predictions at all,
  # we (somewhat arbitrarily) set its
  # precision to 1
  ppv[is.na(ppv)] <- 1
  if (verbose > 2)
    message("+ finished evaluate.get.pr")
  return(list(recall = tpr, precision = ppv))
}

evaluate.calc.aupr <- function(eval,
                               max.tpr = 1,
                               verbose = 0) {
  if (verbose > 2)
    message("+ starting evaluate.calc.aupr")
  pr <- evaluate.get.pr(eval, verbose = verbose)
  idx <- pr$recall <= max.tpr
  if (verbose > 2)
    message("+ finished evaluate.calc.aupr")
  return(evaluate.area.trapez(pr$recall[idx], pr$precision[idx]))
}

create.mlr.learner <- function(method, nrow.data, param.set=NULL,
                               type='BINARY'){
  if (!method %in% c("lasso", "enet", "ridge", "lasso_ll",
                     "ridge_ll", "randomForest", "SVM", "LSVM")){
    stop("Unsupported method!")
  }
  standard.param.set <- list(
    "cost" = c(-2, 3),
    "epsilon" = 1e-08,
    "ntree" = c(100, 1000),
    "mtry" = c(round(sqrt(nrow.data) / 2),
               round(sqrt(nrow.data)),
               round(sqrt(nrow.data) * 2)),
    "alpha" = c(0, 1),
    "class.weights" = c("-1"=5, "1"=1))
  
  if (is.null(param.set)){
    use.param.set <- standard.param.set
  } else {
    use.param.set <- param.set
    for (i in names(standard.param.set)){
      if (!i %in% names(param.set)){
        use.param.set[[i]] <- standard.param.set[[i]]
      }
    }
  }
  
  # lasso/enet/ridge
  if (method %in% c('lasso', 'enet', 'ridge')){
    if (type == 'BINARY'){
      learner <- lrn("classif.cv_glmnet")
      learner$predict_type <- 'prob'
    } else if (type == 'CONTINUOUS'){
      learner <- lrn("regr.cv_glmnet")
    }
    if (method == 'lasso'){
      learner$param_set$values$alpha <- 1
      if ('alpha' %in% names(param.set)){
        warning("Parameter 'alpha' will be ignored and set to 1!")
      }
    } else if (method == 'ridge'){
      learner$param_set$values$alpha <- 0
      if ('alpha' %in% names(param.set)){
        warning("Parameter 'alpha' will be ignored and set to 0!")
      }
    } else if (method == 'enet'){
      if (length(use.param.set$alpha)==1){
        learner$param_set$values$alpha <- use.param.set$alpha
      } else {
        learner$param_set$values$alpha <- to_tune(
          lower=use.param.set$alpha[1],
          upper=use.param.set$alpha[2])
      }
    }
    use.param.set$alpha <- NULL
  } else if (method=='randomForest'){
    if (type == 'BINARY'){
      learner <- lrn("classif.ranger", importance='permutation')
      learner$predict_type <- 'prob'
    } else if (type == 'CONTINUOUS'){
      learner <- lrn("regr.ranger", importance='impurity')
    }
    # number of trees
    if (length(use.param.set$ntree) == 1){
      learner$param_set$values$num.trees <- use.param.set$ntree
    } else if (length(use.param.set$ntree) == 2) {
      learner$param_set$values$num.trees <- to_tune(
        lower=use.param.set$ntree[1], upper=use.param.set$ntree[2])
    } else {
      learner$param_set$values$num.trees <- to_tune(use.param.set$ntree)
    }
    # mtry
    if (length(use.param.set$mtry) == 1){
      learner$param_set$values$mtry <- use.param.set$mtry
    } else if (length(use.param.set$mtry) == 2){
      learner$param_set$values$mtry <- to_tune(
        lower=use.param.set$mtry[1], upper=use.param.set$mtry[2])
    } else {
      learner$param_set$values$mtry <- to_tune(use.param.set$mtry)
    }
    use.param.set$ntree <- NULL
    use.param.set$mtry <- NULL
    if (!'alpha' %in% names(param.set)){
      use.param.set$alpha <- NULL
    }
    if (!'class.weights' %in% names(param.set)){
      use.param.set$class.weights <- NULL
    }
    
  } else if (method %in% c('lasso_ll', 'ridge_ll')){
    if (type == 'BINARY'){
      if (method == 'lasso_ll'){
        type <- 6
      } else {
        type <- 0
      }
      mlr_learners$add("classif.liblinear", LearnerClassifLiblineaR)
      learner <- lrn('classif.liblinear', type=type, 
                     wi=use.param.set$class.weights,
                     epsilon=use.param.set$epsilon)
      learner$predict_type <- 'prob'
    } else if (type == 'CONTINUOUS'){
      stop("Methods not usable for regression tasks!")
    }
    if (length(use.param.set$cost) == 1){
      learner$param_set$values$cost <- use.param.set$cost
    } else if (length(use.param.set$cost) == 2){
      learner$param_set$values$cost <- to_tune(p_dbl(
        lower=use.param.set$cost[1],
        upper=use.param.set$cost[2], trafo=.f_exp))
    } else {
      learner$param_set$values$cost <- to_tune(use.param.set$cost)
    }
    use.param.set$class.weights <- NULL
    use.param.set$epsilon <- NULL
    use.param.set$cost <- NULL
  } else if (method == 'SVM'){
    if (type == 'BINARY'){
      learner <- lrn('classif.svm', kernel='radial')
      learner$predict_type <- 'prob'
      learner$param_set$values$type <- 'C-classification'
      learner$param_set$values$cost <- to_tune(p_dbl(
        lower=-5,
        upper=16, trafo=.f_exp))
      use.param.set$cost <- NULL
      use.param.set$class.weights <- NULL
    } else if (type == 'CONTINUOUS'){
      stop("Methods not usable for regression tasks!")
    }
  } else if (method == 'LSVM'){
    if (type == 'BINARY'){
      learner <- lrn('classif.svm', kernel='linear')
      learner$predict_type <- 'prob'
      learner$param_set$values$type <- 'C-classification'
      learner$param_set$values$cost <- to_tune(p_dbl(
        lower=-5,
        upper=16, trafo=.f_exp))
      use.param.set$cost <- NULL
      use.param.set$class.weights <- NULL
    } else if (type == 'CONTINUOUS'){
      stop("Methods not usable for regression tasks!")
    }
  }
  
  # try to set additional parameters, i hope that mlr catches errors here
  param.settable <- learner$param_set$ids()
  for (x in intersect(names(use.param.set), param.settable)){
    if (length(use.param.set[[x]])==1){
      learner$param_set$values[[x]] <- use.param.set[[x]]
    } else {
      learner$param_set$values[[x]] <- to_tune(use.param.set[[x]])
    }
  }
  
  return(learner)
  
}

perform.feature.selection <- function(data, train.label, param.fs, verbose){
  
  stopifnot(all(c('method', 'no_features', 'direction') %in%names(param.fs)))
  
  # test method.fs
  if (is.factor(train.label)){
    allowed.methods <- c('Wilcoxon', 'AUC', 'gFC')
    if (!param.fs$method %in% allowed.methods) {
      stop('Unrecognised feature selection method. ',
           'Must be one of those: {"',
           paste(allowed.methods, collapse = '", "'), '"}')
    }
  } else if (is.numeric(train.label)){
    allowed.methods <- c('spearman', 'pearson', 'MI')
    if (!param.fs$method %in% allowed.methods) {
      stop('Unrecognised feature selection method. ',
           'Must be one of those: {"',
           paste(allowed.methods, collapse = '", "'), '"}')
    }
  }
  
  # assert the threshold
  stopifnot(param.fs$no_features > 10)
  stopifnot(param.fs$no_features < ncol(data))
  
  if (param.fs$method == 'Wilcoxon') {
    assoc <- vapply(data,
                    FUN=function(x, label){
                      d <- data.frame(x=x, y=label);
                      t <- wilcox.test(x~y, data=d)
                      return(t$p.val)
                    }, FUN.VALUE=double(1),
                    label=train.label)
    assoc <- sort(assoc)
  } else if (param.fs$method == 'AUC') {
    assoc <- vapply(data,
                    FUN=get.single.feat.AUC,
                    FUN.VALUE = double(1),
                    label=train.label,
                    pos=max(levels(train.label)),
                    neg=min(levels(train.label)))
    if (param.fs$direction == 'absolute'){
      assoc[assoc < 0.5] <- 1 - assoc[assoc < 0.5]
    } else if (param.fs$direction == 'negative'){
      assoc <- 1 - assoc
    }
    assoc <- assoc[assoc > 0.5]
    assoc <- sort(-assoc)
  } else if (param.fs$method == 'gFC') {
    assoc <- vapply(data,
                    FUN=get.quantile.FC,
                    FUN.VALUE = double(1),
                    label=train.label,
                    pos=max(levels(train.label)),
                    neg=min(levels(train.label)))
    if (param.fs$direction == 'absolute'){
      assoc <- abs(assoc)
    } else if (param.fs$direction == 'negative'){
      assoc <- -assoc
    }
    assoc <- assoc[assoc > 0]
    assoc <- sort(-assoc)
  } else if (param.fs$method %in% c('spearman', 'pearson')){
    assoc <- vapply(data, FUN=cor, FUN.VALUE = double(1), y=train.label,
                    method=param.fs$method)
    if (param.fs$direction == 'absolute'){
      assoc <- abs(assoc)
    } else if (param.fs$direction == 'negative'){
      assoc <- -assoc
    }
    assoc <- assoc[assoc > 0]
    assoc <- sort(-assoc)
  } else if (param.fs$method == 'MI'){
    assoc <- vapply(data, FUN=function(x){
      mutinformation(discretize(x, disc='equalwidth'),
                     discretize(train.label, disc='equalwidth'))
    }, FUN.VALUE = double(1))
    assoc <- sort(-assoc)
  }
  
  data <- data[,names(assoc)[seq_len(param.fs$no_features)]]
  
  stopifnot(ncol(data) > 0)
  if (verbose > 2) {
    msg <- paste0('++ retaining ', ncol(data),
                  ' features after selection based on ',
                  param.fs$method, '; target number of features ',
                  param.fs$no_features)
    message(msg)
  }
  return(data)
}


measureAUPRC <- function(probs, truth, negative, positive) {
  pr <- pr.curve(scores.class0 = probs[which(truth == positive)],
                 scores.class1 = probs[which(truth == negative)])
  return(pr$auc.integral)
}


get.single.feat.AUC <- function(x, label, pos, neg) {
  x.p <- x[label == pos]
  x.n <- x[label == neg]
  temp.auc <- roc(cases=x.p, controls=x.n, direction='<')$auc
  return(temp.auc)
}


get.quantile.FC <- function(x, label, pos, neg){
  x.p <- x[label == pos]
  x.n <- x[label == neg]
  q.p <- quantile(x.p, probs=seq(.1, .9, length.out=9))
  q.n <- quantile(x.n, probs=seq(.1, .9, length.out=9))
  return(sum(q.p - q.n)/length(q.p))
}

.f_exp <- function(x){10^x}

get.best.glmnet.lambda <- function(model, measure, min.nonzero, task){
  idx <- which(model$model$nzero >= min.nonzero)
  new.model <- model$clone(deep = TRUE)
  perf <- vapply(idx, FUN=function(x){
    new.model$param_set$values$s <- new.model$model$lambda[x]
    pred <- new.model$predict(task)
    pred$score(msr(measure))
  }, FUN.VALUE = double(1))
  if (msr(measure)$minimize){
    f.idx <- which(perf==min(perf))[1]
  } else {
    f.idx <- which(perf==max(perf))[1]
  }
  new.model$param_set$values$s <- new.model$model$lambda[f.idx]
  return(new.model)
}

.rarefy <- function (otu.tab, ss) {  
  set.seed(ss)
  depth = min(rowSums(otu.tab))  
  otu.tab <- as.matrix(otu.tab)  
  ind <- (rowSums(otu.tab) < depth)  
  sam.discard <- rownames(otu.tab)[ind]  
  otu.tab <- otu.tab[!ind, ]  
  rarefy <- function(x, depth) {    
    y <- sample(rep(1:length(x), x), depth)    
    y.tab <- table(y)    
    z <- numeric(length(x))    
    z[as.numeric(names(y.tab))] <- y.tab    
    z  
  }  
  otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))  
  rownames(otu.tab.rff) <- rownames(otu.tab)  
  colnames(otu.tab.rff) <- colnames(otu.tab)  
  return(otu.tab.rff)
}
