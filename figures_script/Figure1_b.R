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
  load("./CRC_Real/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  tax.names <- colnames(data.rel[[1]]$Y)

  # signature effect size
  data.loc <- "./Simulation/By_effsz/"
  facts <- c("sig1.5", "sig2", "sig2.5", "sig3", "sig3.5", "sig4")
  
  PRC_all <- data.frame()
  for(tag in facts){
    load(paste0(data.loc, tag, "/AUPRC/PRC.Rdata"))
    tmp.prc <- PRC %>% group_by(method, linetype) %>% summarize(AUPRC = mean(AUPRC))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label[PRC_all$x.label == "sig1.5"] <- "0.5"
  PRC_all$x.label[PRC_all$x.label == "sig2"] <- "1"
  PRC_all$x.label[PRC_all$x.label == "sig2.5"] <- "1.5"
  PRC_all$x.label[PRC_all$x.label == "sig3"] <- "2"
  PRC_all$x.label[PRC_all$x.label == "sig3.5"] <- "2.5"
  PRC_all$x.label[PRC_all$x.label == "sig4"] <- "3"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$method <- factor(PRC_all$method, levels = unique(PRC_all$method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  
  # Customizing the output
  pdf("./figures/figure1_b_sg_effze.pdf",         # File name
      width = 4.09, height = 3.57, # Width and height in inches
      bg = "white")   
  
  ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(method, linetype), linetype = linetype))  +
    xlab("") + xlab(TeX('$\\Delta$')) + ylim(0.1,0.8) +
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
  data.loc <- "./Simulation/By_taxnum/"
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
  pdf("./figures/figure1_b_sg_sparsity.pdf",         # File name
      width = 4.09, height = 3.57, # Width and height in inches
      bg = "white")   
  
  ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(method, linetype), linetype = linetype, shape = linetype))  +
    xlab("") + xlab("d") + ylim(0.4,0.8) +
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
  data.loc <- "./Simulation/By_pos_prop/"
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
  pdf("./figures/figure1_b_sg_effdir.pdf",         # File name
      width = 4.09, height = 3.57, # Width and height in inches
      bg = "white")   
  
  ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(method, linetype), linetype = linetype))  +
    ylim(0, 1) + xlab("pos%") + ylim(0.4,0.7) +
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
  data.loc <- "./Simulation/By_seqdepth/"
  facts <- c("seq1", "seq1_25", "seq1_50","seq1_75", "seq2")
  PRC_all <- data.frame()
  for(tag in facts){
    load(paste0(data.loc, tag, "/AUPRC/PRC.Rdata"))
    tmp.prc <- PRC %>% group_by(method, linetype) %>% summarize(AUPRC = mean(AUPRC))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label[PRC_all$x.label == "seq1"] <- "1"
  PRC_all$x.label[PRC_all$x.label == "seq1_25"] <- "1.25"
  PRC_all$x.label[PRC_all$x.label == "seq1_50"] <- "1.5"
  PRC_all$x.label[PRC_all$x.label == "seq1_75"] <- "1.75"
  PRC_all$x.label[PRC_all$x.label == "seq2"] <- "2"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$method <- factor(PRC_all$method, levels = unique(PRC_all$method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  # Customizing the output
  pdf("./figures/figure1_b_cc_seqdepth.pdf",         # File name
      width = 4.09, height = 3.57, # Width and height in inches
      bg = "white")   
  
  ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(method, linetype), linetype = linetype))  +
    xlab(TeX('$\\mu$')) + ylim(0.4,0.7) +
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
  
  # signature prevalence
  data.loc <- "./Simulation/By_abd_prop/"
  facts <- c("p20", "p30", "p40", "p50", "p60", "p70","p80")
  PRC_all <- data.frame()
  for(tag in facts){
    load(paste0(data.loc, tag, "/AUPRC/PRC.Rdata"))
    tmp.prc <- PRC %>% group_by(method, linetype) %>% summarize(AUPRC = mean(AUPRC))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label[PRC_all$x.label == "p20"] <- "20"
  PRC_all$x.label[PRC_all$x.label == "p30"] <- "30"
  PRC_all$x.label[PRC_all$x.label == "p40"] <- "40"
  PRC_all$x.label[PRC_all$x.label == "p50"] <- "50"
  PRC_all$x.label[PRC_all$x.label == "p60"] <- "60"
  PRC_all$x.label[PRC_all$x.label == "p70"] <- "70"
  PRC_all$x.label[PRC_all$x.label == "p80"] <- "80"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$method <- factor(PRC_all$method, levels = unique(PRC_all$method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  
  # Customizing the output
  pdf("./figures/figure1_b_sg_prev.pdf",         # File name
      width = 4.09, height = 3.57, # Width and height in inches
      bg = "white")   
  
  ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(method, linetype), linetype = linetype))  +
    xlab(TeX('$\\alpha$')) + ylim(0.5,0.7) +
    geom_line(linewidth = 0.7, aes(color=method)) + geom_point(size = 2, aes(color = method), pch = 18) + scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c("Original", "Batch-corrected"), values=c("solid", "dashed")) +
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