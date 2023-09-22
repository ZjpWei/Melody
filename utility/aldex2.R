library("tidyverse")
library("phyloseq")
library("ALDEx2")

run_aldex2 <- function(otu.tab=feature.table, meta = meta.data, formula){
  countdata <- otu.tab
  mm <- model.matrix(formula(paste0("~", formula)), meta)
  aldex.fit <- aldex.clr(countdata, mm, denom="all")
  x.e <- aldex.glm(aldex.fit)
  glm.effect <- aldex.glm.effect(aldex.fit)
  
  ### adjust p-val by BH, aldex.glm only provide holm q-val.
  p.fdr <- apply(x.e[,grepl("pval", colnames(x.e)) & !grepl("pval.holm", colnames(x.e))], 
                 2, p.adjust, method = "BH")
  colnames(p.fdr) <- paste0(colnames(p.fdr), ".BH")
  x.e <- cbind(x.e, p.fdr)
  return(list(glmfit = x.e, effect = glm.effect))
}
