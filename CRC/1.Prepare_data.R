# =============================================== #
#           CRC 1. CRC data processing            #
# =============================================== #


  # General ----
  rm(list = ls())
  data.loc <- "./CRC/"

  # Prepare data ----
  load("./Data/CRC_data/data/count.Rdata")
  load("./Data/CRC_data/data/meta.Rdata")
  load("./Data/CRC_data/data/tax.Rdata")
  
  meta <- as.data.frame(meta)
  study <- c("AT-CRC", "CN-CRC", "DE-CRC", "FR-CRC", "US-CRC")
  rownames(meta) <- meta$Sample_ID
  meta$Group <- factor(meta$Group, level = c("CTR", "CRC"))
  meta$Study <- factor(meta$Study, levels = study)

  ## CRC1: Austria (109)
  ## CRC2: China (127)
  ## CRC3: German (120)
  ## CRC4: France (114)
  ## CRC5: United States (104)

  ## Remove 2 samples which's sequence depth < 2000
  ## Sample id: "CCIS12370844ST-4-0" and "MMRS51728985ST-27-0-0"
  sample.id.kp <- names(which(rowSums(count) >= 2000))
  meta <- meta[sample.id.kp,]
  count <- count[sample.id.kp,]

  ## CRC1: Austria (CTR: 63, CRC: 46)
  ## CRC2: China (CTR: 54, CRC: 73)
  ## CRC3: German (CTR: 60, CRC: 60)
  ## CRC4: France (CTR: 61, CRC: 52)
  ## CRC5: United States (CTR: 52, CRC: 51)

  ## Batch-corrected by MMUPHin
  batch_adj <- MMUPHin::adjust_batch(feature_abd = t(count),
                                     batch = "Study",
                                     covariates = "Group",
                                     data = meta)
  
  count.batch <- t(batch_adj$feature_abd_adj)
  
  # CRC data K849 ----
  data.rel <- list()
  data.rel.batch <- list()
  for(l in 1:length(study)){
    data.rel[[l]] <- list(Y = count[meta$Sample_ID[meta$Study == study[l]],],
                          X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"))

    data.rel.batch[[l]] <- list(Y = count.batch[meta$Sample_ID[meta$Study == study[l]],],
                                X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"))
  }
  save(data.rel, file = paste0(data.loc, "Processed_data/data.org.K849.Rdata"))
  save(data.rel.batch, file = paste0(data.loc,  "Processed_data/data.bat.K849.Rdata"))

  # CRC data K401 ----
  ## Remove taxa non-zero proportion less than 20%, 401 taxa remain after filtering
  filters <- colMeans(count > 0) >= 0.2
  data.rel <- list()
  data.rel.batch <- list()
  for(l in 1:length(study)){
    data.rel[[l]] <- list(Y = count[meta$Sample_ID[meta$Study == study[l]],filters],
                          X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"))
    data.rel.batch[[l]] <- list(Y = count.batch[meta$Sample_ID[meta$Study == study[l]],filters],
                                X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"))
  }
  save(data.rel, file = paste0(data.loc, "Processed_data/data.org.K401.Rdata"))
  save(data.rel.batch, file = paste0(data.loc,  "Processed_data/data.bat.K401.Rdata"))

  # CRC data K267 (order "186802 Clostridiales") ----
  filters.order <- (tax$order == "186802 Clostridiales") & filters
  data.rel <- list()
  data.rel.batch <- list()
  for(l in 1:length(study)){
    data.rel[[l]] <- list(Y = count[meta$Sample_ID[meta$Study == study[l]],filters.order],
                          X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"))

    data.rel.batch[[l]] <- list(Y = count.batch[meta$Sample_ID[meta$Study == study[l]],filters.order],
                                X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"))
  }
  save(data.rel, file = paste0(data.loc, "Processed_data/data.org.order.Rdata"))
  save(data.rel.batch, file = paste0(data.loc,  "Processed_data/data.bat.order.Rdata"))

