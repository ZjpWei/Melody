  ##############################################
#                                                #
#           Real data processing script          #
#                                                #
  ##############################################
  
  rm(list = ls())
  # set document location
  data.loc <- "./"
  if(!dir.exists(paste0(data.loc, "CRC_Real/CRC_all_K849/prepare_data"))){
    dir.create(paste0(data.loc, "CRC_Real/CRC_all_K849/prepare_data"))
  }
  if(!dir.exists(paste0(data.loc, "CRC_Real/CRC_all_K401/prepare_data"))){
    dir.create(paste0(data.loc, "CRC_Real/CRC_all_K401/prepare_data"))
  }
  if(!dir.exists(paste0(data.loc, "CRC_Real/CRC_all_order/prepare_data"))){
    dir.create(paste0(data.loc, "CRC_Real/CRC_all_order/prepare_data"))
  }
  if(!dir.exists(paste0(data.loc, "CRC_Real/CRC_loso_K849/prepare_data"))){
    dir.create(paste0(data.loc, "CRC_Real/CRC_loso_K849/prepare_data"))
  }
  if(!dir.exists(paste0(data.loc, "CRC_Real/CRC_loso_K401/prepare_data"))){
    dir.create(paste0(data.loc, "CRC_Real/CRC_loso_K401/prepare_data"))
  }
  if(!dir.exists(paste0(data.loc, "CRC_Real/CRC_loso_order/prepare_data"))){
    dir.create(paste0(data.loc, "CRC_Real/CRC_loso_order/prepare_data"))
  }
  if(!dir.exists(paste0(data.loc, "Prediction/LOSO/prepare_data"))){
    dir.create(paste0(data.loc, "Prediction/LOSO/prepare_data"))
  }
  if(!dir.exists(paste0(data.loc, "Prediction/Random_split/prepare_data"))){
    dir.create(paste0(data.loc, "Prediction/Random_split/prepare_data"))
  }
  if(!dir.exists(paste0(data.loc, "figures"))){
    dir.create(paste0(data.loc, "figures"))
  }
  
  # data prepare 
  load("./data/count.Rdata")
  load("./data/meta.Rdata")
  load("./data/tax.Rdata")
  meta <- as.data.frame(meta)
  study <- c("AT-CRC", "CN-CRC", "DE-CRC", "FR-CRC", "US-CRC")
  rownames(meta) <- meta$Sample_ID
  meta$Group <- factor(meta$Group, level = c("CTR", "CRC"))
  meta$Study <- factor(meta$Study, levels = study)
  
  # CRC1: Austria (109)
  # CRC2: China (127)
  # CRC3: German (120)
  # CRC4: France (114)
  # CRC5: United States (104)
  
  # Remove 2 samples which's sequence depth < 2000
  # Sample id: "CCIS12370844ST-4-0" and "MMRS51728985ST-27-0-0"
  sample.id.kp <- names(which(rowSums(count) >= 2000))
  meta <- meta[sample.id.kp,]
  count <- count[sample.id.kp,]
  
  # CRC1: Austria (CTR: 63, CRC: 46)
  # CRC2: China (CTR: 54, CRC: 73)
  # CRC3: German (CTR: 60, CRC: 60)
  # CRC4: France (CTR: 61, CRC: 52)
  # CRC5: United States (CTR: 52, CRC: 51)
  
  # Batch-corrected data
  batch_adj <- MMUPHin::adjust_batch(feature_abd = t(count),
                                     batch = "Study",
                                     covariates = "Group",
                                     data = meta)
  count.batch <- t(batch_adj$feature_abd_adj)
  
  # Save all data for K849
  data.rel <- list()
  data.rel.batch <- list()
  for(l in 1:length(study)){
    data.rel[[l]] <- list(Y = count[meta$Sample_ID[meta$Study == study[l]],],
                          X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"),
                          block = meta$block[meta$Study == study[l]])
    
    data.rel.batch[[l]] <- list(Y = count.batch[meta$Sample_ID[meta$Study == study[l]],],
                                X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"),
                                block = meta$block[meta$Study == study[l]])
  }
  save(data.rel, file = paste0(data.loc, "CRC_Real/CRC_all_K849/prepare_data/data.rel.all.Rdata"))
  save(data.rel.batch, file = paste0(data.loc,  "CRC_Real/CRC_all_K849/prepare_data/data.rel.batch.all.Rdata"))
  
  # Save LOSO data for K849
  tmp.data <- data.rel
  tmp.data.batch <- data.rel.batch
  for(l in 1:length(study)){
    data.rel <- tmp.data[-l]
    data.rel.batch <- tmp.data.batch[-l]
    save(data.rel, file = paste0(data.loc, "CRC_Real/CRC_loso_K849/prepare_data/data.rel.", 
                                 as.character(l), ".Rdata"))
    save(data.rel.batch, file = paste0(data.loc,  "CRC_Real/CRC_loso_K849/prepare_data/data.rel.batch.",
                                       as.character(l), ".Rdata"))
  }
  
  # Save all data for K401
  # remove taxa non-zero proportion less than 20%, 401 taxa remain after filtering
  filters <- colMeans(count > 0) >= 0.2
  data.rel <- list()
  data.rel.batch <- list()
  data.rarify <- list()
  for(l in 1:length(study)){
    data.rel[[l]] <- list(Y = count[meta$Sample_ID[meta$Study == study[l]],filters],
                          X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"),
                          block = meta$block[meta$Study == study[l]])
    data.rel.batch[[l]] <- list(Y = count.batch[meta$Sample_ID[meta$Study == study[l]],filters],
                                X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"),
                                block = meta$block[meta$Study == study[l]])
  }
  save(data.rel, file = paste0(data.loc, "CRC_Real/CRC_all_K401/prepare_data/data.rel.all.Rdata"))
  save(data.rel.batch, file = paste0(data.loc,  "CRC_Real/CRC_all_K401/prepare_data/data.rel.batch.all.Rdata"))
  
  save(data.rel, file = paste0(data.loc, "Prediction/CRC_ltcvo_K401/prepare_data/data.rel.all.Rdata"))
  save(data.rel.batch, file = paste0(data.loc,  "Prediction/CRC_ltcvo_K401/prepare_data/data.rel.batch.all.Rdata"))
  
  save(data.rel, file = paste0(data.loc, "Prediction/CRC_locvo_K401/prepare_data/data.rel.all.Rdata"))
  save(data.rel.batch, file = paste0(data.loc,  "Prediction/CRC_locvo_K401/prepare_data/data.rel.batch.all.Rdata"))
  
  # Save LOSO data for K401
  tmp.data <- data.rel
  tmp.data.batch <- data.rel.batch
  for(l in 1:length(study)){
    data.rel <- tmp.data[-l]
    data.rel.batch <- tmp.data.batch[-l]
    save(data.rel, file = paste0(data.loc, "CRC_Real/CRC_loso_K401/prepare_data/data.rel.", 
                                 as.character(l), ".Rdata"))
    save(data.rel.batch, file = paste0(data.loc,  "CRC_Real/CRC_loso_K401/prepare_data/data.rel.batch.",
                                       as.character(l), ".Rdata"))
  }
  
  # Save data for order "186802 Clostridiales"
  filters.order <- (tax$order == "186802 Clostridiales") & filters
  data.rel <- list()
  data.rel.batch <- list()
  for(l in 1:length(study)){
    data.rel[[l]] <- list(Y = count[meta$Sample_ID[meta$Study == study[l]],filters.order],
                          X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"),
                          block = meta$block[meta$Study == study[l]])
    
    data.rel.batch[[l]] <- list(Y = count.batch[meta$Sample_ID[meta$Study == study[l]],filters.order],
                                X = as.numeric(meta$Group[meta$Study == study[l]] == "CRC"),
                                block = meta$block[meta$Study == study[l]])
  }
  save(data.rel, file = paste0(data.loc, "CRC_Real/CRC_all_order/prepare_data/data.rel.all.Rdata"))
  save(data.rel.batch, file = paste0(data.loc,  "CRC_Real/CRC_all_order/prepare_data/data.rel.batch.all.Rdata"))
  
  # Save LOSO data for order "186802 Clostridiales"
  tmp.data <- data.rel
  tmp.data.batch <- data.rel.batch
  for(l in 1:length(study)){
    data.rel <- tmp.data[-l]
    data.rel.batch <- tmp.data.batch[-l]
    save(data.rel, file = paste0(data.loc, "CRC_Real/CRC_loso_order/prepare_data/data.rel.", 
                                 as.character(l), ".Rdata"))
    save(data.rel.batch, file = paste0(data.loc,  "CRC_Real/CRC_loso_order/prepare_data/data.rel.batch.",
                                       as.character(l), ".Rdata"))
  }
