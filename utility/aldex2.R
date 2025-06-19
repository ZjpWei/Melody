library("tidyverse")
library("phyloseq")
library("ALDEx2")

run_aldex2 <- function(otu.tab=feature.table, meta = meta.data, formula, gamma = NULL){
  countdata <- otu.tab
  mm <- model.matrix(formula(paste0("~", formula)), meta)


  if(is.null(gamma)){
    message("++ Building aldex.clr object without gamma.")
  }else{
    message(paste0("++ Building aldex.clr object with gamma = ", gamma, "."))
  }
  aldex.fit <- aldex.clr(countdata, mm, denom="all", gamma = gamma)

  message("++ Calculate glm test statistics.")
  x.e <- aldex.glm(aldex.fit)

  message("++ Calculate effect sizes and difference between all constrast in glm model.")
  glm.effect <- aldex.glm.effect(aldex.fit)

  ### adjust p-val by BH, aldex.glm only provide holm q-val.
  p.fdr <- apply(x.e[,grepl("pval", colnames(x.e)) & !grepl("pval.holm", colnames(x.e))],
                 2, p.adjust, method = "BH")
  colnames(p.fdr) <- paste0(colnames(p.fdr), ".BH")
  x.e <- cbind(x.e, p.fdr)
  return(list(glmfit = x.e, effect = glm.effect))
}
