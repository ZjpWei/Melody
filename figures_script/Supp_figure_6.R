  ##############################################
#                                                #
#             Supplementary figure 6             #
#                                                #
  ##############################################

  library(phyloseq)
  library(ggplot2)
  
  rm(list = ls())
  set.seed(2023)
  load("./data/meta.Rdata")
  load("./CRC_Real/CRC_all_K849/prepare_data/data.rel.all.Rdata")
  source("./utility/GDM_utility.R")
  
  # Prepare plootting
  title <- "CRC Simulated data"
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
    colnames(dmData) <- colnames(data.rel[[l]]$Y)
    
    Simulate.count.rel = NULL
    for(ll in 1:nrow(dmData)){
      Simulate.count.rel <- cbind(Simulate.count.rel, rmultinom(1, 1e5, dmData[ll,]))
    }
    feature.table <- rbind(feature.table, t(Simulate.count.rel)/rowSums(t(Simulate.count.rel)))
    outcome <- c(outcome, data.rel[[l]]$X)
    Study <- c(Study, rep(paste0("Study",as.character(l)), length(data.rel[[l]]$X)))
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
  pdf("./figures/Supp_6.pdf",         # File name
      width =  6.80, height = 4.66, # Width and height in inches
      bg = "white")   
  
  ggplot(PCoA, aes(x=PCo1, y=PCo2, color = Study)) + geom_point() +
    theme(text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14),
          legend.position = "right") + guides(color=guide_legend(title="Simulated study"))
  
  
  # Closing the graphical device
  dev.off() 
  