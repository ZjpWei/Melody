  ##############################################
#                                                #
#         Sensitivity: Melody models with        #
#           different reference taxon            #
#                                                #
  ##############################################
  
  library("miMeta")
  rm(list = ls())
  load("./CRC_Real/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  L <- length(data.rel)
  
  # Different reference taxa for Melody models
  # ref1: "Coprococcus catus [ref_mOTU_v2_4874]"
  # ref2: "unknown Lachnospiraceae [meta_mOTU_v2_6937]"
  # ref3: "unknown Butyricicoccus [meta_mOTU_v2_5437]"
  # ref4: "Dorea formicigenerans [ref_mOTU_v2_0973]"  
  # ref5: "unknown Clostridiales [meta_mOTU_v2_6832]"
  
  Study <- c()
  Group <- c()
  Block <- c()
  rel.abd <- c()
  for(l in 1:length(data.rel)){
    Study <- c(Study, rep(as.character(l), nrow(data.rel[[l]]$Y)))
    Group <- c(Group, as.character(data.rel[[l]]$X))
    Block <- c(Block, factor(data.rel[[l]]$block))
    rel.abd <- rbind(rel.abd, data.rel[[l]]$Y)
  }
  rel.abd <- t(rel.abd)
  meta <- data.frame(Sample_ID = colnames(rel.abd), Study = Study, Group = Group, Block = Block)
  rownames(meta) <- meta$Sample_ID
  ######################## Generate summary statsitics #########################
  # ref1
  summary.stat.study.ref1 <- melody.get.summary(rel.abd = rel.abd,
                                                sample.data = meta,
                                                sample.id = "Sample_ID",
                                                study = "Study",
                                                disease = "Group",
                                                ref = "Coprococcus catus [ref_mOTU_v2_4874]")
  
  save(summary.stat.study.ref1, file = "./CRC_Real/CRC_all_K401/summary.stats/summary.stats.ref1.Rdata")
  # ref2
  summary.stat.study.ref2 <- melody.get.summary(rel.abd = rel.abd,
                                                sample.data = meta,
                                                sample.id = "Sample_ID",
                                                study = "Study",
                                                disease = "Group",
                                                ref = "unknown Lachnospiraceae [meta_mOTU_v2_6937]")
  
  save(summary.stat.study.ref2, file = "./CRC_Real/CRC_all_K401/summary.stats/summary.stats.ref2.Rdata")
  # ref3
  summary.stat.study.ref3 <- melody.get.summary(rel.abd = rel.abd,
                                                sample.data = meta,
                                                sample.id = "Sample_ID",
                                                study = "Study",
                                                disease = "Group",
                                                ref = "unknown Butyricicoccus [meta_mOTU_v2_5437]")
  
  save(summary.stat.study.ref3, file = "./CRC_Real/CRC_all_K401/summary.stats/summary.stats.ref3.Rdata")
  #ref4
  summary.stat.study.ref4 <- melody.get.summary(rel.abd = rel.abd,
                                                sample.data = meta,
                                                sample.id = "Sample_ID",
                                                study = "Study",
                                                disease = "Group",
                                                ref = "Dorea formicigenerans [ref_mOTU_v2_0973]")
  
  save(summary.stat.study.ref4, file = "./CRC_Real/CRC_all_K401/summary.stats/summary.stats.ref4.Rdata")
  # ref5
  summary.stat.study.ref5 <- melody.get.summary(rel.abd = rel.abd,
                                                sample.data = meta,
                                                sample.id = "Sample_ID",
                                                study = "Study",
                                                disease = "Group",
                                                ref = "unknown Clostridiales [meta_mOTU_v2_6832]")
  
  save(summary.stat.study.ref5, file = "./CRC_Real/CRC_all_K401/summary.stats/summary.stats.ref5.Rdata")
  
  summary.stat.study.mixref <- melody.get.summary(rel.abd = rel.abd,
                                                  sample.data = meta,
                                                  sample.id = "Sample_ID",
                                                  study = "Study",
                                                  disease = "Group",
                                                  ref = c("Coprococcus catus [ref_mOTU_v2_4874]",
                                                          "unknown Lachnospiraceae [meta_mOTU_v2_6937]",
                                                          "unknown Butyricicoccus [meta_mOTU_v2_5437]",       
                                                          "Dorea formicigenerans [ref_mOTU_v2_0973]",
                                                          "unknown Clostridiales [meta_mOTU_v2_6832]"))
  
  save(summary.stat.study.mixref, file = "./CRC_Real/CRC_all_K401/summary.stats/summary.stat.study.mixref.Rdata")
  ################################ Run Melody models #############################
  # Melody.model.ref1
  Melody.model.ref1 <- melody.meta.summary(Melody = summary.stat.study.ref1, tune.type = "HBIC")
  
  save(Melody.model.ref1, file = "./CRC_Real/CRC_all_K401/sensitivity/Melody.model.ref1.Rdata")
  
  # Melody.model.ref2
  Melody.model.ref2 <- melody.meta.summary(Melody = summary.stat.study.ref2, tune.type = "HBIC")
  
  save(Melody.model.ref2, file = "./CRC_Real/CRC_all_K401/sensitivity/Melody.model.ref2.Rdata")
  
  # Melody.model.ref3
  Melody.model.ref3 <- melody.meta.summary(Melody = summary.stat.study.ref3, tune.type = "HBIC")
  
  save(Melody.model.ref3, file = "./CRC_Real/CRC_all_K401/sensitivity/Melody.model.ref3.Rdata")
  
  # Melody.model.ref4
  Melody.model.ref4 <- melody.meta.summary(Melody = summary.stat.study.ref4, tune.type = "HBIC")
  
  save(Melody.model.ref4, file = "./CRC_Real/CRC_all_K401/sensitivity/Melody.model.ref4.Rdata")
  
  # Melody.model.ref5
  Melody.model.ref5 <- melody.meta.summary(Melody = summary.stat.study.ref5, tune.type = "HBIC")
  
  save(Melody.model.ref5, file = "./CRC_Real/CRC_all_K401/sensitivity/Melody.model.ref5.Rdata")
  
  # Melody.model.mixref
  Melody.model.mixref <- melody.meta.summary(Melody = summary.stat.study.mixref, tune.type = "HBIC")
  
  save(Melody.model.mixref, file = "./CRC_Real/CRC_all_K401/sensitivity/Melody.model.mixref.Rdata")
