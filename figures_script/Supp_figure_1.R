# =============================================== #
#             Supplementary figure 1              #
# =============================================== #

  # Packages ----
  library("phyloseq")
  library("ggplot2")

  # General ----
  rm(list = ls())
  
  ## Original ----
  load("./CRC_meta_analysis/data/meta.Rdata")
  load("./CRC_meta_analysis/CRC_all_K849/prepare_data/data.rel.all.Rdata")
  title <- "CRC Real"
  studymeta <- unique(meta$Study)
  feature.table <- NULL
  outcome <- NULL
  Study <- NULL
  K <- 849
  for(l in 1:length(data.rel)){
    feature.table <- rbind(feature.table, data.rel[[l]]$Y )
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

  ### Generate figures ----
  pdf("./figures/Supp_1_PCoA_raw.pdf", width = 6.14, height = 4.66, bg = "white")

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

  dev.off()

  ## Batch-corrected ----
  rm(list = ls())
  load("./CRC_meta_analysis/data/meta.Rdata")
  load("./CRC_meta_analysis/CRC_all_K849/prepare_data/data.rel.batch.all.Rdata")
  title <- "batch corrected"

  studymeta <- unique(meta$Study)
  feature.table <- NULL
  outcome <- NULL
  Study <- NULL
  K <- 849
  for(l in 1:length(data.rel.batch)){
    feature.table <- rbind(feature.table, data.rel.batch[[l]]$Y)
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

  ### Generate figures ----
  pdf("./figures/Supp_1_PCoA_batch.pdf", width = 6.14, height = 4.66, bg = "white")

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

  dev.off()
