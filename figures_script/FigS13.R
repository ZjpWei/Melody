# =============================================== #
#             Supplementary figure 13             #
# =============================================== #
  
  # Packages ----
  library('ggplot2')
  library("tidyverse")
  library("latex2exp")
  library("ggtext")
  
  # Main ----
  rm(list = ls())
  load("./CRC/Processed_data/data.org.K401.Rdata")
  tax.names <- colnames(data.rel[[1]]$Y)
  
  # Scenario: n = 100 ~ 180 ----
  ## Signature effect size ----
  data.loc <- "./Simulation/large/Independent/"
  facts <- c("1", "1.5", "2", "2.5", "3")
  PRC_all <- NULL
  for(tag in facts){
    PRC_tmp <- NULL
    for(s in 1:100){
      load(paste0(data.loc, "AUPRC_Ka40_Pos0.7_effsz", tag, "_mu0_", as.character(s), ".Rdata"))
      PRC_tmp <- rbind(PRC_tmp, PRC)
    }
    tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180 <br> Signature effect size"
  
  p11 <- PRC_all %>% filter(method != "Melody-Pooled") %>% 
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab(TeX('$\\Delta$', bold = TRUE)) + ylim(0.39,0.83) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) + 
    theme_bw() +
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
    labs(linetype="Data") +
    facet_wrap(facets = vars(scenarios)) + guides(colour = guide_legend(nrow = 1))
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180<br> Signature effect size"
  
  ## Signature sparsity ----
  facts <- c("20", "30", "40", "50","60", "70", "80")
  PRC_all <- NULL
  for(tag in facts){
    PRC_tmp <- NULL
    for(s in 1:100){
      load(paste0(data.loc, "AUPRC_Ka", tag, "_Pos0.7_effsz2_mu0_", as.character(s), ".Rdata"))
      PRC_tmp <- rbind(PRC_tmp, PRC)
    }
    tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180 <br> Signature number"
  
  p12 <-  PRC_all %>% filter(method != "Melody-Pooled") %>% 
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab(TeX("d", italic = TRUE)) + ylim(0.39,0.83) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) +
    theme_bw() +
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
    facet_wrap(facets = vars(scenarios)) + guides(colour = guide_legend(nrow = 1))
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180<br> Signature number"
  
  ## Signature effect direction ----
  facts <- c("50","60", "70", "80", "90","100")
  PRC_all <- NULL
  for(tag in facts){
    PRC_tmp <- NULL
    for(s in 1:100){
      load(paste0(data.loc, "AUPRC_Ka40_Pos", as.numeric(tag)/100, "_effsz2_mu0_", as.character(s), ".Rdata"))
      PRC_tmp <- rbind(PRC_tmp, PRC)
    }
    tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180 <br> Signature effect direction"
  
  p13 <- PRC_all %>% filter(method != "Melody-Pooled") %>%
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab(TeX("pos%", italic = TRUE)) + ylim(0.39,0.83) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) +
    theme_bw() +
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
    facet_wrap(facets = vars(scenarios)) + guides(colour = guide_legend(nrow = 1))
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180<br> Signature effect direction"
  
  ## Case/control sequence depth unevenness ----
  facts <- c("0", "0.25", "0.50","0.75", "1")
  PRC_all <- NULL
  for(tag in facts){
    PRC_tmp <- NULL
    for(s in 1:100){
      load(paste0(data.loc, "AUPRC_Ka40_Pos0.7_effsz2_mu", tag, "_", as.character(s), ".Rdata"))
      PRC_tmp <- rbind(PRC_tmp, PRC)
    }
    tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180 <br> Sequencing depth unevenness"
  
  p14 <- PRC_all %>% filter(method != "Melody-Pooled") %>%
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab(TeX('u', italic = TRUE)) + ylim(0.39,0.83)+ 
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2, aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) + 
    theme_bw() +
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
    facet_wrap(facets = vars(scenarios)) + guides(colour = guide_legend(nrow = 1))
  PRC_all$scenarios <- "<i>n</i> = 100 ~ 180<br> Sequencing depth unevenness"
  
  # Scenario: n = 20 ~ 60 ----
  ## Signature effect size ----
  data.loc <- "./Simulation/small/Independent/"
  facts <- c("4", "4.5", "5", "5.5", "6")
  PRC_all <- NULL
  for(tag in facts){
    PRC_tmp <- NULL
    for(s in 1:100){
      load(paste0(data.loc, "AUPRC_Ka40_Pos0.7_effsz", tag, "_mu0_", as.character(s), ".Rdata"))
      PRC_tmp <- rbind(PRC_tmp, PRC)
    }
    tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60 <br> Signature effect size"
  
  p21 <- PRC_all %>% filter(method != "Melody-Pooled") %>% 
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab(TeX('$\\Delta$', bold = TRUE)) + ylim(0.39,0.83) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) + 
    theme_bw() +
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
    labs(linetype="Data") +
    facet_wrap(facets = vars(scenarios)) + guides(colour = guide_legend(nrow = 1))
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60<br> Signature effect size"
  
  ## Signature sparsity ----
  facts <- c("20", "30", "40", "50","60", "70", "80")
  PRC_all <- NULL
  for(tag in facts){
    PRC_tmp <- NULL
    for(s in 1:100){
      load(paste0(data.loc, "AUPRC_Ka", tag, "_Pos0.7_effsz2_mu0_", as.character(s), ".Rdata"))
      PRC_tmp <- rbind(PRC_tmp, PRC)
    }
    tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60 <br> Signature number"
  
  p22 <-  PRC_all %>% filter(method != "Melody-Pooled") %>% 
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab(TeX("d", italic = TRUE)) + ylim(0.39,0.83) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) +
    theme_bw() +
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
    facet_wrap(facets = vars(scenarios)) + guides(colour = guide_legend(nrow = 1))
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60<br> Signature number"
  
  ## Signature effect direction ----
  facts <- c("50","60", "70", "80", "90","100")
  PRC_all <- NULL
  for(tag in facts){
    PRC_tmp <- NULL
    for(s in 1:100){
      load(paste0(data.loc, "AUPRC_Ka40_Pos", as.numeric(tag)/100, "_effsz2_mu0_", as.character(s), ".Rdata"))
      PRC_tmp <- rbind(PRC_tmp, PRC)
    }
    tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60 <br> Signature effect direction"
  
  p23 <- PRC_all %>% filter(method != "Melody-Pooled") %>%
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab(TeX("pos%", italic = TRUE)) + ylim(0.39,0.83) +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) +
    theme_bw() +
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
    facet_wrap(facets = vars(scenarios)) + guides(colour = guide_legend(nrow = 1))
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60<br> Signature effect direction"
  
  ## Case/control sequence depth unevenness ----
  facts <- c("0", "0.25", "0.50","0.75", "1")
  PRC_all <- NULL
  for(tag in facts){
    PRC_tmp <- NULL
    for(s in 1:100){
      load(paste0(data.loc, "AUPRC_Ka40_Pos0.7_effsz2_mu", tag, "_", as.character(s), ".Rdata"))
      PRC_tmp <- rbind(PRC_tmp, PRC)
    }
    tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
    tmp.prc$x.label <- tag
    PRC_all <- rbind(PRC_all, tmp.prc)
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60 <br> Sequencing depth unevenness"
  
  p24 <- PRC_all %>% filter(method != "Melody-Pooled") %>%
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab(TeX('u', italic = TRUE)) + ylim(0.39,0.83)+ 
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2, aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) + 
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="bottom",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16)) +
    facet_wrap(facets = vars(scenarios)) + guides(colour = guide_legend(nrow = 1))
  PRC_all$scenarios <- "<i>n</i> = 20 ~ 60<br> Sequencing depth unevenness"
  
  pdf("./figures/FigS13.pdf", width = 15, height = 7.5, bg = "white")
  
  ggpubr::annotate_figure(
    ggpubr::ggarrange(p11, p12, p13, p14, p21, p22, p23, p24, nrow = 2, ncol = 4,
                      common.legend = TRUE, legend = "bottom", label.y = "AUPRC"),
    left = grid::textGrob("AUPRC", rot = 90, hjust = -0.5, vjust = 0.5, 
                          gp = grid::gpar(cex = 1.6))
  )
  
  dev.off()
