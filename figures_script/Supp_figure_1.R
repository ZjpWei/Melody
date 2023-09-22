  ##############################################
#                                                #
#             Supplementary figure 1             #
#                                                #
  ##############################################
  
  library("phyloseq")
  library("ggplot2")
  
  rm(list = ls())
  # General
  set.seed(2023)
  source("./utility/GDM_utility.R")
  load("./data/meta.Rdata")
  load("./CRC_Real/CRC_all_K849/prepare_data/data.rel.all.Rdata")
  title <- "CRC Real"
  studymeta <- unique(meta$Study)
  feature.table <- NULL
  outcome <- NULL
  Study <- NULL
  K <- 849
  for(l in 1:length(data.rel)){
    load(paste0("./GDM_fit/GDM.", as.character(l),".Rdata"))
    X <- matrix(1, nrow = nrow(data.rel[[l]]$Y), ncol = 1)
    dmData <- r.GDM(X = X, mod = mod.gdm)
    Y.tmp <- data.rel[[l]]$Y
    idx <- order(colMeans(Y.tmp/rowSums(Y.tmp)), decreasing = TRUE)
    idx.rev <- rank(-colMeans(Y.tmp/rowSums(Y.tmp)))
    stopifnot(sum(colMeans(Y.tmp/rowSums(Y.tmp)) != 0) == ncol(dmData))
    dmData <- cbind(dmData, matrix(0, nrow = nrow(data.rel[[l]]$Y), ncol = K - ncol(dmData)))[,idx.rev]
    colnames(dmData) <- colnames( data.rel[[l]]$Y)
    
    feature.table <- rbind(feature.table, dmData)
    outcome <- c(outcome, data.rel[[l]]$X)
    Study <- c(Study, rep(paste0("CRC",as.character(l)), length(data.rel[[l]]$X)))
  }
  rownames(feature.table) <- paste0("sa", as.character(1:length(outcome)))
  names(outcome) <- rownames(feature.table)
  meta.data = data.frame(labels = outcome)
  rownames(feature.table) <- rownames(meta.data)
  
  feature.table = t(feature.table)
  OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
  META = sample_data(meta.data)
  PHYSEQ = phyloseq(OTU, META)
  
  # Calculate ordination
  iMDS <- ordinate(physeq = PHYSEQ, "PCoA")
  
  PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"], PCo2 = iMDS$vectors[,"Axis.2"], Study = Study,
                     outcome = as.character(outcome), variable = title)
  
  # Customizing the output
  pdf("./figures/Supp_1_PCoA_raw.pdf",         # File name
      width = 6.14, height = 4.66, # Width and height in inches
      bg = "white")   
  
  ggplot(PCoA, aes(x=PCo1, y=PCo2, color=Study)) + geom_point() +
    theme(text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.key = element_rect(fill = "white"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.position = "right")
  
  # Closing the graphical device
  dev.off() 
  ##################################. batch. #####################################
  
  rm(list = ls())
  set.seed(2023)
  load("~/Documents/Melody/data/meta.Rdata")
  load("~/Documents/Melody/CRC_K849/prepare_data/data.rel.batch.Rdata")
  
  title <- "batch corrected"
  source("./utility/GDM_utility.R")
  
  studymeta <- unique(meta$Study)
  feature.table <- NULL
  outcome <- NULL
  Study <- NULL
  K <- 849
  for(l in 1:length(data.rel.batch)){
    load(paste0("./GDM_fit/GDM.batch.", as.character(l),".Rdata"))
    X <- matrix(1, nrow = nrow(data.rel.batch[[l]]$Y), ncol = 1)
    dmData <- r.GDM(X = X, mod = mod.gdm)
    Y.tmp <-  data.rel.batch[[l]]$Y
    idx <- order(colMeans(Y.tmp/rowSums(Y.tmp)), decreasing = TRUE)
    idx.rev <- rank(-colMeans(Y.tmp/rowSums(Y.tmp)))
    stopifnot(sum(colMeans(Y.tmp/rowSums(Y.tmp)) != 0) == ncol(dmData))
    dmData <- cbind(dmData, matrix(0, nrow = nrow(data.rel.batch[[l]]$Y), ncol = K - ncol(dmData)))[,idx.rev]
    colnames(dmData) <- colnames( data.rel.batch[[l]]$Y)
    
    feature.table <- rbind(feature.table, dmData)
    outcome <- c(outcome, data.rel.batch[[l]]$X)
    Study <- c(Study, rep(paste0("CRC",as.character(l)), length(data.rel.batch[[l]]$X)))
  }
  
  rownames(feature.table) <- paste0("sa", as.character(1:length(outcome)))
  names(outcome) <- rownames(feature.table)
  meta.data = data.frame(labels = outcome)
  rownames(feature.table) <- rownames(meta.data)
  
  feature.table = t(feature.table)
  OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
  META = sample_data(meta.data)
  PHYSEQ = phyloseq(OTU, META)
  
  # Calculate ordination
  iMDS <- ordinate(physeq = PHYSEQ, "PCoA")
  
  PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"], PCo2 = iMDS$vectors[,"Axis.2"], Study = Study,
                     outcome = as.character(outcome), variable = title)
  
  # Customizing the output
  pdf("./figures/Supp_1_PCoA_batch.pdf",         # File name
      width = 6.14, height = 4.66, # Width and height in inches
      bg = "white")   
  
  ggplot(PCoA, aes(x=PCo1, y=PCo2, color=Study)) + geom_point() +
    theme(text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.key = element_rect(fill = "white"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.position = "right") 
  
  # Closing the graphical device
  dev.off() 