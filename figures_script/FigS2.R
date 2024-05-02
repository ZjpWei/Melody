# =============================================== #
#             Supplementary figure 2              #
# =============================================== #

  # Packages ----
  library("phyloseq")
  library("ggplot2")

  # General ----
  rm(list = ls())

  ## Original ----
  load("./Data/CRC_data/data/meta.Rdata")
  load("./CRC/Processed_data/data.org.K849.Rdata")
  title <- "CRC Real"
  studymeta <- unique(meta$Study)
  feature.table <- NULL
  outcome <- NULL
  Study <- NULL
  K <- 849
  for(l in 1:length(data.rel)){
    feature.table <- rbind(feature.table, data.rel[[l]]$Y/rowSums(data.rel[[l]]$Y) )
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
  iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "bray")

  PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"], PCo2 = iMDS$vectors[,"Axis.2"], Study = Study,
                     outcome = as.character(outcome), variable = title)

  ### Generate figures ----
  pdf("./figures/FigS2_raw.pdf", width = 6.14, height = 4.66, bg = "white")

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

  ## Simulation data ----
  rm(list = ls())
  source("./utility/GDM_utility.R")
  load("./Data/CRC_data/data/meta.Rdata")
  load("./CRC/Processed_data/data.org.K401.Rdata")
  title <- "Simulation data"

  set.seed(2023)
  data.null <- list()
  K <- 401
  n.sample <- c(100, 120, 140, 160, 180)
  Y.study <- list()
  for(l in 1:5){
    load(paste0("./Simulation/GDM_fit/GDM.", as.character(l), ".Rdata"))
    X <- matrix(1, nrow = n.sample[l], ncol = 1)
    dmData <- r.GDM(X = X, mod = mod.gdm)
    Y.study[[l]] <- data.rel[[l]]$Y
    idx <- order(colMeans(Y.study[[l]]/rowSums(Y.study[[l]])), decreasing = TRUE)
    idx.rev <- rank(-colMeans(Y.study[[l]]/rowSums(Y.study[[l]])))
    stopifnot(sum(colMeans(Y.study[[l]]/rowSums(Y.study[[l]])) != 0) == ncol(dmData))
    dmData <- cbind(dmData, matrix(0, nrow = n.sample[l], ncol = K - ncol(dmData)))[,idx.rev]
    colnames(dmData) <- colnames(data.rel[[l]]$Y)
    data.null[[l]] <- list(Y = dmData, X = rep(0, nrow(dmData)))
  }

  data.rel <- list()
  for(l in 1:5){
    c.name <- colnames(Y.study[[l]])
    n = n.sample[l]
    # generate absolute abundant data
    # shuffle subjects (we always assign the first half subjects to cases)
    case.idx.1 = 1:round(n/2)
    Simulate.count.1 = data.null[[l]]$Y
    Prob.abs.1 = Simulate.count.1/rowSums(Simulate.count.1)

    # generate relative abundant data
    Simulate.depth.1 <- sample(x = rowSums(Y.study[[l]]), size = n, replace = TRUE)
    Simulate.count.rel.1 = NULL
    for(ll in 1:nrow(Simulate.count.1)){
      Simulate.count.rel.1 <- cbind(Simulate.count.rel.1, rmultinom(1, Simulate.depth.1[ll], Prob.abs.1[ll,]))
    }
    Simulate.count.rel.1 <- t(Simulate.count.rel.1)
    # generate taxa and sample names
    data.rel[[l]] <- list(Y=Simulate.count.rel.1,X = rep(0, nrow(Simulate.count.rel.1)))
    sample.nm <- paste0("Cohort.",as.character(l),".s",as.character(1:n))
    colnames(data.rel[[l]]$Y) <- c.name
    rownames(data.rel[[l]]$Y) <- sample.nm
  }

  studymeta <- unique(meta$Study)
  feature.table <- NULL
  outcome <- NULL
  Study <- NULL
  for(l in 1:length(data.rel)){
    feature.table <- rbind(feature.table, data.rel[[l]]$Y/rowSums(data.rel[[l]]$Y))
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
  iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "bray")

  PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"], PCo2 = iMDS$vectors[,"Axis.2"], Study = Study,
                     outcome = as.character(outcome), variable = title)

  ### Generate figures ----
  pdf("./figures/FigS2_MelodySIM.pdf", width = 6.14, height = 4.66, bg = "white")

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

  # ## Batch-corrected ----
  # rm(list = ls())
  # load("./CRC_meta_analysis/data/meta.Rdata")
  # load("./CRC_meta_analysis/CRC_all_K849/prepare_data/data.rel.batch.all.Rdata")
  # title <- "batch corrected"
  #
  # studymeta <- unique(meta$Study)
  # feature.table <- NULL
  # outcome <- NULL
  # Study <- NULL
  # K <- 849
  # for(l in 1:length(data.rel.batch)){
  #   feature.table <- rbind(feature.table, data.rel.batch[[l]]$Y/rowSums(data.rel.batch[[l]]$Y))
  #   outcome <- c(outcome, data.rel.batch[[l]]$X)
  #   Study <- c(Study, rep(paste0("CRC",as.character(l)), length(data.rel.batch[[l]]$X)))
  # }
  #
  # rownames(feature.table) <- paste0("sa", as.character(1:length(outcome)))
  # names(outcome) <- rownames(feature.table)
  # meta.data = data.frame(labels = outcome)
  # rownames(feature.table) <- rownames(meta.data)
  #
  # feature.table = t(feature.table)
  # OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
  # META = sample_data(meta.data)
  # PHYSEQ = phyloseq(OTU, META)
  #
  # # Calculate ordination
  # iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "bray")
  #
  # PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"], PCo2 = iMDS$vectors[,"Axis.2"], Study = Study,
  #                    outcome = as.character(outcome), variable = title)
  #
  # ### Generate figures ----
  # pdf("./figures/Supp_figure1_PCoA_batch.pdf", width = 6.14, height = 4.66, bg = "white")
  #
  # ggplot(PCoA, aes(x=PCo1, y=PCo2, color=Study)) + geom_point() +
  #   theme(text = element_text(size = 16),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         panel.background = element_rect(fill = 'white'),
  #         panel.border = element_rect(colour = "black", fill=NA),
  #         legend.key = element_rect(fill = "white"),
  #         legend.title = element_text(size = 14),
  #         legend.text = element_text(size = 13),
  #         legend.position = "right")
  #
  # dev.off()
  #
  # ## ANCOMBC2 Simulation data ----
  # rm(list = ls())
  # load("./CRC_meta_analysis/data/meta.Rdata")
  # load("./CRC_meta_analysis/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  # source("./utility/GDM_utility.R")
  # title <- "simulation data"
  #
  # set.seed(2023)
  # data.null <- list()
  # K <- 401
  # n.sample <- c(100, 120, 140, 160, 180)
  # for(l in 1:5){
  #   if(l != 2){
  #     X <- matrix(1, nrow = n.sample[l], ncol = 1)
  #     Y.tmp <- data.rel[[l]]$Y
  #     dmData <- ANCOMBC::sim_plnm(abn_table = Y.tmp, taxa_are_rows = FALSE,
  #                                 prv_cut = 0, n = n.sample[l], lib_mean = 1e8,
  #                                 disp = 0.5)
  #
  #     log_abn_data = log(dmData + 1e-5) ## Weizhou comments: count data + 1e-5
  #     n <- n.sample[l]
  #     d <- nrow(log_abn_data)
  #     rownames(log_abn_data) = paste0("T", seq_len(d))
  #     colnames(log_abn_data) = paste0("S", seq_len(n))
  #
  #     # Generate the sample and feature meta data
  #     # Sampling fractions are set to differ by the variable of interest
  #     smd = data.frame(sample = paste0("S", seq_len(n)),
  #                      samp_frac = log(c(runif(n/2, min = 1e-4, max = 1e-3),
  #                                        runif(n/2, min = 1e-3, max = 1e-2))),
  #                      cont_cov = rnorm(n),
  #                      bin_cov = as.factor(rep(seq_len(2), each = n/2)))
  #
  #
  #     fmd = data.frame(taxon = paste0("T", seq_len(d)),
  #                      seq_eff = log(runif(d, min = 0.1, max = 1)))
  #
  #     # Add effect sizes of covariates to the true abundances
  #     # smd_dmy = model.matrix(~ 0 + cont_cov + bin_cov, data = smd)
  #     # log_abn_data = log_abn_data + outer(fmd$lfc_cont, smd_dmy[, "cont_cov"] )
  #     # log_abn_data = log_abn_data + outer(fmd$lfc_bin, smd_dmy[, "bin_cov2"])
  #
  #     # Add sample- and taxon-specific biases
  #     log_otu_data = t(t(log_abn_data) + smd$samp_frac)
  #     log_otu_data = log_otu_data + fmd$seq_eff
  #
  #     ## Weizhou comments: no round
  #     otu_data = round(exp(log_otu_data))
  #     dmData <- t(as.matrix(otu_data))
  #     colnames(dmData) <- colnames(data.rel[[l]]$Y)
  #     dmData <- dmData[rowSums(dmData) > 0,]
  #     data.null[[l]] <- list(Y = dmData, X = rep(0, nrow(dmData)))
  #   }
  # }
  #
  # studymeta <- unique(meta$Study)
  # feature.table <- NULL
  # outcome <- NULL
  # Study <- NULL
  # for(l in 1:length(data.null)){
  #   if(l != 2){
  #     feature.table <- rbind(feature.table, data.null[[l]]$Y/rowSums(data.null[[l]]$Y))
  #     outcome <- c(outcome, data.null[[l]]$X)
  #     Study <- c(Study, rep(paste0("CRC",as.character(l)), length(data.null[[l]]$X)))
  #   }
  # }
  # rownames(feature.table) <- paste0("sa", as.character(1:length(outcome)))
  # names(outcome) <- rownames(feature.table)
  # meta.data = data.frame(labels = outcome)
  # rownames(feature.table) <- rownames(meta.data)
  #
  # feature.table = t(feature.table)
  # OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
  # META = sample_data(meta.data)
  # PHYSEQ = phyloseq(OTU, META)
  #
  # # Calculate ordination
  # iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "bray")
  #
  # PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"], PCo2 = iMDS$vectors[,"Axis.2"], Study = Study,
  #                    outcome = as.character(outcome), variable = title)
  #
  # ### Generate figures ----
  # pdf("./figures/Supp_PCoA_ANCOMSIM.pdf", width = 6.14, height = 4.66, bg = "white")
  #
  # ggplot(PCoA, aes(x=PCo1, y=PCo2, color=Study)) + geom_point() +
  #   theme(text = element_text(size = 16),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         panel.background = element_rect(fill = 'white'),
  #         panel.border = element_rect(colour = "black", fill=NA),
  #         legend.key = element_rect(fill = "white"),
  #         legend.title = element_text(size = 14),
  #         legend.text = element_text(size = 13),
  #         legend.position = "right")
  #
  # dev.off()
