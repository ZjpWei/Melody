# =============================================== #
#     CRC 4.Sensitivity: Melody models with       #
#            different reference taxon            #
# =============================================== #

  # Packages ----
  library("miMeta")

  # General ----
  rm(list = ls())
  load("./CRC/Processed_data/data.org.K401.Rdata")
  L <- length(data.rel)

  ## Different reference taxa for Melody models
  ## ref1: "Coprococcus catus [ref_mOTU_v2_4874]"
  ## ref2: "unknown Lachnospiraceae [meta_mOTU_v2_6937]"
  ## ref3: "unknown Butyricicoccus [meta_mOTU_v2_5437]"
  ## ref4: "Dorea formicigenerans [ref_mOTU_v2_0973]"
  ## ref5: "unknown Clostridiales [meta_mOTU_v2_7130]"

  rel.abd <- list()
  covariate.interest <- list()
  for(l in 1:L){
    rel.abd[[paste0("S",l)]] <-  data.rel[[l]]$Y
    covariate.interest[[paste0("S",l)]] <- data.frame("disease" = data.rel[[l]]$X)
  }
  ref.1 <- rep("Coprococcus catus [ref_mOTU_v2_4874]",L)
  ref.2 <- rep("unknown Lachnospiraceae [meta_mOTU_v2_6937]",L)
  ref.3 <- rep("unknown Butyricicoccus [meta_mOTU_v2_5437]",L)
  ref.4 <- rep("Dorea formicigenerans [ref_mOTU_v2_0973]",L)
  ref.5 <- rep("unknown Clostridiales [meta_mOTU_v2_7130]",L)
  ref.mix <- c("unknown Lachnospiraceae [meta_mOTU_v2_6937]",
               "unknown Clostridiales [meta_mOTU_v2_7130]",
               "unknown Butyricicoccus [meta_mOTU_v2_5437]",
               "Coprococcus catus [ref_mOTU_v2_4874]",
               "Dorea formicigenerans [ref_mOTU_v2_0973]")
  
  names(ref.1) <- names(rel.abd)
  names(ref.2) <- names(rel.abd)
  names(ref.3) <- names(rel.abd)
  names(ref.4) <- names(rel.abd)
  names(ref.5) <- names(rel.abd)
  names(ref.mix) <- names(rel.abd)
  
  # Generate summary statistics ----
  ## ref1
  null.obj.1 <- melody.null.model(rel.abd = rel.abd, ref = ref.1)
  
  summary.stat.study.ref1 <- melody.get.summary(null.obj = null.obj.1, covariate.interest = covariate.interest)

  save(summary.stat.study.ref1, file = "./CRC/Sensitivity/summary.stats.ref1.Rdata")
  
  ## ref2
  null.obj.2 <- melody.null.model(rel.abd = rel.abd, ref = ref.2)
  
  summary.stat.study.ref2 <- melody.get.summary(null.obj = null.obj.2, covariate.interest = covariate.interest)

  save(summary.stat.study.ref2, file = "./CRC/Sensitivity/summary.stats.ref2.Rdata")
  
  ## ref3
  null.obj.3 <- melody.null.model(rel.abd = rel.abd, ref = ref.3)
  
  summary.stat.study.ref3 <- melody.get.summary(null.obj = null.obj.3, covariate.interest = covariate.interest)
  
  save(summary.stat.study.ref3, file = "./CRC/Sensitivity/summary.stats.ref3.Rdata")

  ## ref4
  null.obj.4 <- melody.null.model(rel.abd = rel.abd, ref = ref.4)
  
  summary.stat.study.ref4 <- melody.get.summary(null.obj = null.obj.4, covariate.interest = covariate.interest)
  
  save(summary.stat.study.ref4, file = "./CRC/Sensitivity/summary.stats.ref4.Rdata")

  ## ref5
  null.obj.5 <- melody.null.model(rel.abd = rel.abd, ref = ref.5)
  
  summary.stat.study.ref5 <- melody.get.summary(null.obj = null.obj.5, covariate.interest = covariate.interest)
  
  save(summary.stat.study.ref5, file = "./CRC/Sensitivity/summary.stats.ref5.Rdata")
  
  ## ref.mix
  null.obj.mix <- melody.null.model(rel.abd = rel.abd, ref = ref.mix)
  
  summary.stat.study.mixref <- melody.get.summary(null.obj = null.obj.mix, covariate.interest = covariate.interest)
  
  save(summary.stat.study.mixref, file = "./CRC/Sensitivity/summary.stat.study.mixref.Rdata")
  
  # Run Melody models ----
  ## Melody model with reference 1
  Melody.model.ref1 <- melody.meta.summary(summary.stats = summary.stat.study.ref1)

  save(Melody.model.ref1, file = "./CRC/Sensitivity/Melody.model.ref1.Rdata")

  ## Melody model with reference 2
  Melody.model.ref2 <- melody.meta.summary(summary.stats = summary.stat.study.ref2)

  save(Melody.model.ref2, file = "./CRC/Sensitivity/Melody.model.ref2.Rdata")

  ## Melody model with reference 3
  Melody.model.ref3 <- melody.meta.summary(summary.stats = summary.stat.study.ref3)

  save(Melody.model.ref3, file = "./CRC/Sensitivity/Melody.model.ref3.Rdata")

  ## Melody model with reference 4
  Melody.model.ref4 <- melody.meta.summary(summary.stats = summary.stat.study.ref4)

  save(Melody.model.ref4, file = "./CRC/Sensitivity/Melody.model.ref4.Rdata")

  ## Melody model with reference 5
  Melody.model.ref5 <- melody.meta.summary(summary.stats = summary.stat.study.ref5)

  save(Melody.model.ref5, file = "./CRC/Sensitivity/Melody.model.ref5.Rdata")

  ## Melody model with mixed references
  Melody.model.mixref <- melody.meta.summary(summary.stats = summary.stat.study.mixref)

  save(Melody.model.mixref, file = "./CRC/Sensitivity/Melody.model.mixref.Rdata")
