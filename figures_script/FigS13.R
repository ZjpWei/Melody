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
  data.loc <- "./Simulation/large/Independent/"
  facts <- c("4", "20", "40", "80", "160")

  PRC_all <- NULL
  for(pos.lst in c(0.6, 0.8, 1)){
    for(u.lst in c(0, 0.5, 1)){
      for(tag in facts){
        PRC_tmp <- NULL
        for(s in 1:100){
          load(paste0(data.loc, "AUPRC_Ka", tag, "_Pos", pos.lst, "_effsz2_mu", u.lst, "_", as.character(s), ".Rdata"))
          PRC_tmp <- rbind(PRC_tmp, PRC)
        }
        tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, batch_type, data_type) %>%
          summarize(mean_auprc = mean(AUPRC, na.rm = TRUE), error = sd(AUPRC, na.rm = TRUE)) %>%
          dplyr::transmute(AUPRC = mean_auprc, stderr = error)
        tmp.prc$upper <- pmin(tmp.prc$AUPRC + tmp.prc$stderr, 1)
        tmp.prc$lower <- pmax(tmp.prc$AUPRC - tmp.prc$stderr, 0.2)
        tmp.prc$x.label <- tag
        tmp.prc$pos <- paste0("<i>pos</i>% = ", pos.lst * 100, "%")
        tmp.prc$u <- paste0("<i>u</i> = ", u.lst)
        PRC_all <- rbind(PRC_all, tmp.prc)
      }
    }
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$method %in% c("Melody", "Melody-median"),]
  PRC_all$Method <- factor(PRC_all$method, levels = c("Melody", "Melody-median"), ordered = TRUE)
  PRC_all <- PRC_all[PRC_all$batch_type %in% c("Original","MMUPHin"),]
  PRC_all$batch_type[PRC_all$batch_type == "MMUPHin"] <- "Study-effect corrected"
  PRC_all$pos <- factor(PRC_all$pos, levels = paste0("<i>pos</i>% = ",  c(0.6, 0.8, 1) * 100, "%"))

  p.13 <- PRC_all %>% filter(method != "Melody-Pooled") %>%
    ggplot(aes(x=x.label, y=AUPRC, group = interaction(Method, batch_type), linetype = batch_type))  +
    xlab("Number of signatures") + ylab("AUPRC") + ylim(0.1, 1) + ggtitle("Sample size 100-180") +
    geom_line(linewidth = 0.7, aes(color=Method), position = position_dodge(width = 0.2)) +
    geom_point(size = 2,aes(color = Method), position = position_dodge(width = 0.2), pch = 18) +
    scale_color_manual(breaks =  c("Melody", "Melody-median"), values = c("red","pink"))  +
    scale_linetype_manual(breaks = c( "Original", "Study-effect corrected"), values=c("solid", "dashed")) +
    geom_errorbar(aes(ymin = AUPRC - stderr, ymax = AUPRC + stderr, color = Method),width = 0.1,
                  position = position_dodge(width = 0.2)) + theme_bw() +
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

  pdf("./figures/FigS13.pdf", width = 8.50, height = 8.00, bg = "white")

  p.13

  dev.off()
