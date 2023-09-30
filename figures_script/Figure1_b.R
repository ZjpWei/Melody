  ##############################################
#                                                #
#    Figure 1 (b) Simulation AUPRC Comparison    #
#                                                #
  ##############################################

  library('ggplot2')
  library("tidyverse")
  library("latex2exp")

  rm(list = ls())
  # General
  load("./CRC_meta_analysis/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  tax.names <- colnames(data.rel[[1]]$Y)

  # signature effect size
  data.loc <- "./Simulation/Sig_effsz/"
  facts <- c("sig1", "sig1.5", "sig2", "sig2.5", "sig3")

  PRC_all <- data.frame()
  for(tag in facts){
    load(paste0(data.loc, tag, "/AUPRC/PRC.Rdata"))
    tmp.prc <- PRC %>% group_by(method, linetype) %>% summarize(AUPRC = mean(AUPRC))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label[PRC_all$x.label == "sig1"] <- "1"
  PRC_all$x.label[PRC_all$x.label == "sig1.5"] <- "1.5"
  PRC_all$x.label[PRC_all$x.label == "sig2"] <- "2"
  PRC_all$x.label[PRC_all$x.label == "sig2.5"] <- "2.5"
  PRC_all$x.label[PRC_all$x.label == "sig3"] <- "3"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$method <- factor(PRC_all$method, levels = unique(PRC_all$method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]

  # Customizing the output
  pdf("./figures/figure1_b_sig_effsz.pdf",         # File name
      width = 4.09, height = 3.57, # Width and height in inches
      bg = "white")

  ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(method, linetype), linetype = linetype))  +
    xlab("") + xlab(TeX('$\\Delta$')) + ylim(0.3,0.8) +
    geom_line(linewidth = 0.7, aes(color=method)) + geom_point(size = 2,aes(color = method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 14),
          legend.key =  element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "none",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.75, 'cm')) +
    labs(linetype="Pooled data")

  # Closing the graphical device
  dev.off()

  # signature sparsity
  data.loc <- "./Simulation/Sig_number/"
  facts <- c("d20", "d30", "d40", "d50","d60", "d70","d80")

  PRC_all <- data.frame()
  for(tag in facts){
    load(paste0(data.loc, tag, "/AUPRC/PRC.Rdata"))
    tmp.prc <- PRC %>% group_by(method, linetype) %>% summarize(AUPRC = mean(AUPRC))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label[PRC_all$x.label == "d20"] <- "20"
  PRC_all$x.label[PRC_all$x.label == "d30"] <- "30"
  PRC_all$x.label[PRC_all$x.label == "d40"] <- "40"
  PRC_all$x.label[PRC_all$x.label == "d50"] <- "50"
  PRC_all$x.label[PRC_all$x.label == "d60"] <- "60"
  PRC_all$x.label[PRC_all$x.label == "d70"] <- "70"
  PRC_all$x.label[PRC_all$x.label == "d80"] <- "80"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$method <- factor(PRC_all$method, levels = unique(PRC_all$method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  # Customizing the output
  pdf("./figures/figure1_b_sig_number.pdf",         # File name
      width = 4.09, height = 3.57, # Width and height in inches
      bg = "white")

  ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(method, linetype), linetype = linetype, shape = linetype))  +
    xlab("") + xlab("d") + ylim(0.3,0.8) +
    geom_line(linewidth = 0.7, aes(color=method)) + geom_point(size = 2,aes(color = method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) +
    scale_shape_manual(breaks = c( "Original", "Batch-corrected"), values=c(18,18,20)) +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 14))

  # Closing the graphical device
  dev.off()

  # signature effect direction
  data.loc <- "./Simulation/Sig_effdir/"
  facts <- c("p50","p60", "p70", "p80", "p90","p100")
  PRC_all <- data.frame()
  for(tag in facts){
    load(paste0(data.loc, tag, "/AUPRC/PRC.Rdata"))
    tmp.prc <- PRC %>% group_by(method, linetype) %>% summarize(AUPRC = mean(AUPRC))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label[PRC_all$x.label == "p50"] <- "50"
  PRC_all$x.label[PRC_all$x.label == "p60"] <- "60"
  PRC_all$x.label[PRC_all$x.label == "p70"] <- "70"
  PRC_all$x.label[PRC_all$x.label == "p80"] <- "80"
  PRC_all$x.label[PRC_all$x.label == "p90"] <- "90"
  PRC_all$x.label[PRC_all$x.label == "p100"] <- "100"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$method <- factor(PRC_all$method, levels = unique(PRC_all$method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]

  # Customizing the output
  pdf("./figures/figure1_b_sig_effdir.pdf",         # File name
      width = 4.09, height = 3.57, # Width and height in inches
      bg = "white")

  ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(method, linetype), linetype = linetype))  +
    ylim(0, 1) + xlab("pos%") + ylim(0.3,0.8) +
    geom_line(linewidth = 0.7, aes(color=method)) + geom_point(size = 2,aes(color = method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 14))

  # Closing the graphical device
  dev.off()

  # case/control sequence depth unevenness
  data.loc <- "./Simulation/Sig_depth/"
  facts <- c("seq0", "seq0.25", "seq0.50","seq0.75", "seq1")
  PRC_all <- data.frame()
  for(tag in facts){
    load(paste0(data.loc, tag, "/AUPRC/PRC.Rdata"))
    tmp.prc <- PRC %>% group_by(method, linetype) %>% summarize(AUPRC = mean(AUPRC))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label[PRC_all$x.label == "seq0"] <- "0"
  PRC_all$x.label[PRC_all$x.label == "seq0.25"] <- "0.25"
  PRC_all$x.label[PRC_all$x.label == "seq0.50"] <- "0.50"
  PRC_all$x.label[PRC_all$x.label == "seq0.75"] <- "0.75"
  PRC_all$x.label[PRC_all$x.label == "seq1"] <- "1"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$method <- factor(PRC_all$method, levels = unique(PRC_all$method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  # Customizing the output
  pdf("./figures/figure1_b_sig_depth.pdf",         # File name
      width = 4.09, height = 3.57, # Width and height in inches
      bg = "white")

  ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(method, linetype), linetype = linetype))  +
    xlab(TeX('$\\mu$')) + ylim(0.3,0.8) +
    geom_line(linewidth = 0.7, aes(color=method)) + geom_point(size = 2, aes(color = method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 14))

  # Closing the graphical device
  dev.off()
