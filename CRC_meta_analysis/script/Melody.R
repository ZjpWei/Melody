# =============================================== #
#                Meta-analysis                    # 
# =============================================== #

  rm(list = ls())
  
  # Packages ----
  library("miMeta")
  
  # Parameters ----
  seed <- 2023

  # General ----
  cat('Starting model building script\n')
  start.time <- proc.time()[1]

  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("The analysis tag needs to be provided! Exiting...\n")
  }
  s <- args[1]
  tag <- args[2]

  data.loc <- paste0("./CRC_meta_analysis/", tag, "/")
  if(!dir.exists(paste0(data.loc, "Models_original"))){
    dir.create(paste0(data.loc, "Models_original"))
  }
  load(paste0(data.loc, "prepare_data/data.rel.", s, ".Rdata"))
  data.rel.analysis <- data.rel

  # Run model ----
  L <- length(data.rel.analysis)
  Study <- NULL
  Group <- NULL
  rel.abd <- NULL
  for(l in 1:length(data.rel.analysis)){
    Study <- c(Study, rep(as.character(l), nrow(data.rel.analysis[[l]]$Y)))
    Group <- c(Group, as.character(data.rel.analysis[[l]]$X))
    rel.abd <- rbind(rel.abd, data.rel.analysis[[l]]$Y)
  }
  rel.abd <- t(rel.abd)
  meta <- data.frame(Sample_ID = colnames(rel.abd), Study = Study, Group = Group)
  rownames(meta) <- meta$Sample_ID

  # Decide reference
  if(tag == "CRC_all_K401" | tag == "CRC_loso_K401"){
    ref <- "Coprococcus catus [ref_mOTU_v2_4874]"
  }else if(tag == "CRC_all_K849" | tag == "CRC_loso_K849"){
    ref <- "Coprococcus catus [ref_mOTU_v2_4874]"
  }else if(tag == "CRC_all_order" | tag == "CRC_loso_order"){
    ref <- "unknown Lachnospiraceae [meta_mOTU_v2_6937]"
  }
  
  # Generate summary statistics
  summary.stat.study <- melody.get.summary(rel.abd = rel.abd,
                                           sample.data = meta,
                                           sample.id = "Sample_ID",
                                           study = "Study",
                                           disease = "Group",
                                           ref = ref)
  
  # meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study)
  
  save(Melody.model, file = paste0(data.loc, "Models_original/Melody.model.", s, ".Rdata"))
