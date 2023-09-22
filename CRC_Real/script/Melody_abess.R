  ##############################################
#                                                #
#          Meta-Analysis Model Building          #
#                                                #
  ##############################################
  # Packages
  library('coin')
  library('phyloseq')
  library('pracma')
  library('corpcor')
  library('MASS')
  library("miMeta")
  
  # Parameters
  ref.1 <- "Coprococcus catus [ref_mOTU_v2_4874]"
  seed <- 2023
  
  # General 
  cat('Starting model building script\n')
  start.time <- proc.time()[1]
  
  # args = commandArgs(trailingOnly=TRUE)
  # if (length(args)==0) {
  #   stop("The analysis tag needs to be provided! Exiting...\n")
  # }
  s <- "all" #args[1]
  tag <- "CRC_all_order" #args[2]
  
  data.loc <- paste0("./CRC_Real/", tag, "/")
  ###############################################################
  load(paste0(data.loc, "prepare_data/data.rel.", s, ".Rdata"))
  data.rel.analysis <- data.rel
  
  # Run models
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
  
  # generate summary statistics
  summary.stat.study <- melody.get.summary(rel.abd = rel.abd,
                                           sample.data = meta,
                                           sample.id = "Sample_ID",
                                           study = "Study",
                                           disease = "Group",
                                           ref = ref.1,
                                           verbose = TRUE)
  
  # meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study, verbose = TRUE)
  
  save(Melody.model, file = paste0(data.loc, "Models/Melody.model.", s, ".Rdata"))
  