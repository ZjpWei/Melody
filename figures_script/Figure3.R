# =============================================== #
#      Figure 3 Simulation AUPRC Comparison       #
# =============================================== #

  # Packages ----
  library('ggplot2')
  library("tidyverse")
  library("latex2exp")
  library("ggtext")
  
  # General ----
  rm(list = ls())
  load("./CRC_meta_analysis/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  tax.names <- colnames(data.rel[[1]]$Y)
  
  # Scenario: n = 100 ~ 180 ----
  ## Signature effect size ----
  data.loc <- "./Simulation/large/Sig_effsz/"
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
  PRC_all$Method <- factor(PRC_all$Method, levels = unique(PRC_all$Method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$linetype %in% c("Original","MMUPHin"),]
  PRC_all$linetype[PRC_all$linetype == "MMUPHin"] <- "Batch-corrected"
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180 <br> Signature effect size"
  
  p11 <- ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(Method, linetype), linetype = linetype))  +
    xlab(TeX('$\\Delta$', bold = TRUE)) + ylim(0.3,0.81) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) + theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.key =  element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "none",
          legend.box="vertical",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(1.3, 'cm'),
          strip.text = element_markdown(size = 16)) +
    labs(linetype="Pooled data") +
    facet_wrap(facets = vars(scenarios))
  
  ## Signature sparsity ----
  data.loc <- "./Simulation_large/Sig_number/"
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
  PRC_all$Method <- factor(PRC_all$Method, levels = unique(PRC_all$Method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$linetype %in% c("Original","MMUPHin"),]
  PRC_all$linetype[PRC_all$linetype == "MMUPHin"] <- "Batch-corrected"
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180 <br> Signature number"
  
  p12 <- ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(Method, linetype), linetype = linetype, shape = linetype))  +
    xlab(TeX("d", italic = TRUE)) + ylim(0.3,0.81) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) +
    scale_shape_manual(breaks = c( "Original", "Batch-corrected"), values=c(18,18,20)) + theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16)) +
    facet_wrap(facets = vars(scenarios))
  
  ## Signature effect direction ----
  data.loc <- "./Simulation_large/Sig_effdir/"
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
  PRC_all$Method <- factor(PRC_all$Method, levels = unique(PRC_all$Method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$linetype %in% c("Original","MMUPHin"),]
  PRC_all$linetype[PRC_all$linetype == "MMUPHin"] <- "Batch-corrected"
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180 <br> Signature effect direction"
  
  p13 <- ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(Method, linetype), linetype = linetype))  +
    ylim(0, 1) + xlab(TeX("pos%", italic = TRUE)) + ylim(0.3,0.81) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) + theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16)) +
    facet_wrap(facets = vars(scenarios))
  
  ## Case/control sequence depth unevenness ----
  data.loc <- "./Simulation_large/Sig_depth/"
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
  PRC_all$x.label[PRC_all$x.label == "seq0.50"] <- "0.5"
  PRC_all$x.label[PRC_all$x.label == "seq0.75"] <- "0.75"
  PRC_all$x.label[PRC_all$x.label == "seq1"] <- "1"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$Method, levels = unique(PRC_all$Method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$linetype %in% c("Original","MMUPHin"),]
  PRC_all$linetype[PRC_all$linetype == "MMUPHin"] <- "Batch-corrected"
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180 <br> Sequencing depth unevenness"
  
  p14 <- ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(Method, linetype), linetype = linetype))  +
    xlab(TeX('u', italic = TRUE)) + ylim(0.3,0.81) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2, aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) + theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16)) +
    facet_wrap(facets = vars(scenarios))
  
  # Scenario: n = 20 ~ 60 ----
  ## Signature effect size ----
  data.loc <- "./Simulation/small/Sig_effsz/"
  facts <- c("sig4", "sig4.5", "sig5", "sig5.5", "sig6")
  
  PRC_all <- data.frame()
  for(tag in facts){
    load(paste0(data.loc, tag, "/AUPRC/PRC.Rdata"))
    tmp.prc <- PRC %>% group_by(method, linetype) %>% summarize(AUPRC = mean(AUPRC))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label[PRC_all$x.label == "sig4"] <- "4"
  PRC_all$x.label[PRC_all$x.label == "sig4.5"] <- "4.5"
  PRC_all$x.label[PRC_all$x.label == "sig5"] <- "5"
  PRC_all$x.label[PRC_all$x.label == "sig5.5"] <- "5.5"
  PRC_all$x.label[PRC_all$x.label == "sig6"] <- "6"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$Method, levels = unique(PRC_all$Method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$linetype %in% c("Original","MMUPHin"),]
  PRC_all$linetype[PRC_all$linetype == "MMUPHin"] <- "Batch-corrected"
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60 <br> Signature effect size"
  
  p21 <- ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(Method, linetype), linetype = linetype))  +
    xlab(TeX('$\\Delta$', bold = TRUE)) + ylim(0.3,0.8) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) + theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          axis.text = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.key =  element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "none",
          legend.box="vertical",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.75, 'cm'),
          strip.text = element_markdown(size = 16)) +
    labs(linetype="Pooled data") +
    facet_wrap(facets = vars(scenarios))
  
  
  ## Signature sparsity ----
  data.loc <- "./Simulation_samll/Sig_number/"
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
  PRC_all$Method <- factor(PRC_all$Method, levels = unique(PRC_all$Method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$linetype %in% c("Original","MMUPHin"),]
  PRC_all$linetype[PRC_all$linetype == "MMUPHin"] <- "Batch-corrected"
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60 <br> Signature number"
  
  p22 <- ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(Method, linetype), linetype = linetype, shape = linetype))  +
    xlab("<i>d</i>") + ylim(0.3,0.8) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) +
    scale_shape_manual(breaks = c( "Original", "Batch-corrected"), values=c(18,18,20)) + theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          axis.text = element_text(size = 14),
          axis.title.x = element_markdown(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16)) +
    facet_wrap(facets = vars(scenarios))
  
  ## Signature effect direction ----
  data.loc <- "./Simulation_small/Sig_effdir/"
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
  PRC_all$Method <- factor(PRC_all$Method, levels = unique(PRC_all$Method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$linetype %in% c("Original","MMUPHin"),]
  PRC_all$linetype[PRC_all$linetype == "MMUPHin"] <- "Batch-corrected"
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60 <br> Signature effect direction"
  
  p23 <- ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(Method, linetype), linetype = linetype))  +
    xlab("<i>pos%</i>") + ylim(0.3,0.8) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) + theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 14),
          axis.title.x = element_markdown(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16)) +
    facet_wrap(facets = vars(scenarios))
  
  ## Case/control sequence depth unevenness ----
  data.loc <- "./Simulation_small/Sig_depth/"
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
  PRC_all$x.label[PRC_all$x.label == "seq0.50"] <- "0.5"
  PRC_all$x.label[PRC_all$x.label == "seq0.75"] <- "0.75"
  PRC_all$x.label[PRC_all$x.label == "seq1"] <- "1"
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$Method, levels = unique(PRC_all$Method), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody","ALDEx2","ANCOM-BC","BW", "CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$linetype %in% c("Original","MMUPHin"),]
  PRC_all$linetype[PRC_all$linetype == "MMUPHin"] <- "Batch-corrected"
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60 <br> Sequencing depth unevenness"
  
  p24 <- ggplot(PRC_all, aes(x=x.label, y=AUPRC, group = interaction(Method, linetype), linetype = linetype))  +
    xlab("<i>u</i>") + ylim(0.3,0.8) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2, aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(breaks = c( "Original", "Batch-corrected"), values=c("solid", "dashed")) + theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 14),
          axis.title.x = element_markdown(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16)) +
    facet_wrap(facets = vars(scenarios))
  
  pp.all <- ggpubr::ggarrange(p11, p12, p13, p14, p21, p22, p23, p24, nrow = 2, ncol = 4,
                              common.legend = TRUE, legend = "bottom", label.y = "AUPRC")
  
  # Generate figures ----
  pdf("./figures/figure3.pdf", width = 15, height = 7.5, bg = "white")
  
  ggpubr::annotate_figure(pp.all, left = grid::textGrob("AUPRC", rot = 90, vjust = 0.5, hjust = -0.3,
                                                        gp = grid::gpar(cex = 1.5)))
  
  dev.off()
