huge.pr <- function (path, theta, verbose = TRUE, plot = TRUE) {
  gcinfo(verbose = FALSE)
  ROC = list()
  d = length(theta)
  pos.total = sum(theta != 0)
  neg.total = d - pos.total
  if (verbose)
    message("Computing F1 scores, false positive rates and true positive rates....", appendLF=FALSE)
  ROC$prec = rep(0, length(path))
  ROC$rec  = rep(0, length(path))
  ROC$F1 = rep(0, length(path))
  for (r in 1:length(path)) {
    tmp = path[[r]]
    tp.all = (theta != 0) * (tmp != 0)
    ROC$tp[r] <- sum(tp.all != 0)/pos.total
    fp.all = (theta == 0) * (tmp != 0)
    ROC$fp[r] <- sum(fp.all != 0)/neg.total
    fn.all = (theta != 0) * (tmp == 0)
    precision = sum(tp.all)/(sum(tp.all) + sum(fp.all))
    recall = sum(tp.all)/(sum(tp.all) + sum(fn.all))
    
    ROC$prec[r] <- precision
    ROC$rec[r]  <- recall
    ROC$F1[r] = 2 * precision * recall/(precision + recall)
    if (is.na(ROC$F1[r]))
      ROC$F1[r] = 0
  }
  if (verbose)
    message("done.")
  rm(precision, recall, tp.all, fp.all, fn.all, path, theta)
  gc()
  ord.p = order(ROC$prec, ROC$rec, na.last=NA)
  ROC$prec <- ROC$prec[ord.p]
  ROC$rec  <- ROC$rec[ord.p]
  tmp2 = c(min(c(pos.total/d, ROC$prec)), ROC$prec, 1)
  tmp1 = c(1, ROC$rec, 0)
  if (plot) {
    par(mfrow = c(1, 1))
    plot(tmp1, tmp2, type = "b", main = "PR Curve", xlab = "Recall",
         ylab = "Precision", ylim = c(0, 1))
  }
  tmax <- diff(range(tmp2))*diff(range(tmp1))
  ROC$AUC = sum(diff(tmp2) * (tmp1[-1] + tmp1[-length(tmp1)]))/(2*tmax)
  rm(ord.p, tmp1, tmp2)
  gc()
  class(ROC) = "roc"
  return(ROC)
}