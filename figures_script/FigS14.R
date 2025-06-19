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
  data.loc <- "./Simulation/large/Independent_RASim_320/"

  PRC_all <- NULL
  for(pos.lst in c(0.6, 0.8, 1)){
    for(u.lst in c(0, 0.5, 1)){
      PRC_tmp <- NULL
      for(s in 1:100){
        load(paste0(data.loc, "AUPRC_Ka320_Pos", pos.lst, "_effsz2_mu", u.lst, "_", as.character(s), ".Rdata"))
        PRC_tmp <- rbind(PRC_tmp, PRC)
      }
      tmp.prc <- data.frame(PRC_tmp)
      tmp.prc$pos <- paste0("<i>pos</i>% = ", pos.lst * 100, "%")
      tmp.prc$u <- paste0("<i>u</i> = ", u.lst)
      PRC_all <- rbind(PRC_all, tmp.prc)
    }
  }
  PRC_all <- PRC_all[PRC_all$method %in% c("Melody", "ALDEx2-default", "ANCOM-BC2", "MMUPHin","BW","CLR-LASSO"),]
  PRC_all$method[PRC_all$method == "ALDEx2-default"] <- "ALDEx2"
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "MMUPHin", "CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","MMUPHin"),]
  PRC_all$batch_type[PRC_all$batch_type == "MMUPHin"] <- "Study-effect corrected"
  PRC_all$pos <- factor(PRC_all$pos, levels = paste0("<i>pos</i>% = ",  c(0.6, 0.8, 1) * 100, "%"))

  p.large <- PRC_all %>%
    ggplot(aes(x=Method, y=AUPRC, linetype = batch_type))  +
    xlab("Number of signatures") + ylab("AUPRC") + ylim(0.2,1) + ggtitle("Sample size 100-180") +
    geom_boxplot(
      aes(color = Method, linetype = batch_type, group = interaction(Method, batch_type)),
      position = position_dodge(width = 0.8, preserve = "single"),
      linewidth = 0.7,
      width = 0.6
    ) +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red","purple", "orange", "grey", "green" ,"blue"))  +  theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_markdown(hjust = 0.5, size = 20),
          axis.title = element_text(size = 16),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text.y = element_text(size = 14),
          axis.text.x = element_markdown(angle = 30, vjust = 1, hjust=1, size = 14),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          strip.text = element_markdown(size = 16),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(colour = guide_legend(nrow = 1))

  pdf("./figures/FigS14.pdf", width = 8.60, height = 7.84, bg = "white")

  p.large

  dev.off()
