# =============================================== #
#             Supplementary figure 12             #
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
  data.loc <- "./Simulation_v2/large/Independent/"
  facts <- c("20", "30", "40", "50", "60", "70", "80")

  PRC_all <- NULL
  for(pos.lst in c(0.6, 0.8, 1)){
    for(u.lst in c(0, 0.5, 1)){
      for(tag in facts){
        PRC_tmp <- NULL
        for(s in 1:100){
          load(paste0(data.loc, "AUPRC_Ka", tag, "_Pos", pos.lst, "_effsz2_mu", u.lst, "_", as.character(s), ".Rdata"))
          PRC_tmp <- rbind(PRC_tmp, PRC)
        }
        tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
        tmp.prc$x.label <- tag
        tmp.prc$pos <- paste0("<i>pos</i>% = ", pos.lst * 100, "%")
        tmp.prc$u <- paste0("<i>u</i> = ", u.lst)
        PRC_all <- rbind(PRC_all, tmp.prc)
      }
    }
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$pos <- factor(PRC_all$pos, levels = paste0("<i>pos</i>% = ",  c(0.6, 0.8, 1) * 100, "%"))

  p.large <- PRC_all %>% filter(method != "Melody-Pooled") %>%
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab("Number of signatures") + ylim(0.2,0.9) + ggtitle("Sample size 100-180") +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_markdown(hjust = 0.5, size = 20),
          axis.title = element_text(size = 16),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(colour = guide_legend(nrow = 1))

  # Scenario: n = 20 ~ 60 ----
  data.loc <- "./Simulation_v2/small/Independent/"
  facts <- c("20", "30", "40", "50", "60", "70", "80")

  PRC_all <- NULL
  for(pos.lst in c(0.6, 0.8, 1)){
    for(u.lst in c(0, 0.5, 1)){
      for(tag in facts){
        PRC_tmp <- NULL
        for(s in 1:100){
          if(file.exists(paste0(data.loc, "AUPRC_Ka", tag, "_Pos", pos.lst, "_effsz5_mu", u.lst, "_", as.character(s), ".Rdata"))){
            load(paste0(data.loc, "AUPRC_Ka", tag, "_Pos", pos.lst, "_effsz5_mu", u.lst, "_", as.character(s), ".Rdata"))
            PRC_tmp <- rbind(PRC_tmp, PRC)
          }
        }
        tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>% summarize(AUPRC = mean(AUPRC, na.rm = TRUE))
        tmp.prc$x.label <- tag
        tmp.prc$pos <- paste0("<i>pos</i>% = ", pos.lst * 100, "%")
        tmp.prc$u <- paste0("<i>u</i> = ", u.lst)
        PRC_all <- rbind(PRC_all, tmp.prc)
      }
    }
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-Pooled", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$Method %in% c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","ConQuR"),]
  PRC_all$batch_type[PRC_all$batch_type == "ConQuR"] <- "Study-effect corrected"
  PRC_all$pos <- factor(PRC_all$pos, levels = paste0("<i>pos</i>% = ",  c(0.6, 0.8, 1) * 100, "%"))

  p.small <- PRC_all %>% filter(method != "Melody-Pooled") %>%
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab("Number of signatures") + ylim(0.2,0.8) + ggtitle("Sample size 20-60") +
    geom_line(linewidth = 0.7, aes(color=Method)) + geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_markdown(hjust = 0.5, size = 20),
          axis.title = element_text(size = 16),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(colour = guide_legend(nrow = 1))

  pdf("./figures/FigS12.pdf", width = 13.83, height = 9, bg = "white")

  ggpubr::annotate_figure(
    ggpubr::ggarrange(p.large, p.small, nrow = 1, ncol = 2,
                      common.legend = TRUE, legend = "bottom", label.y = "AUPRC"))

  dev.off()
