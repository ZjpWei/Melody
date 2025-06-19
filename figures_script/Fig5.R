# =============================================== #
#          Figure 5: Prediction AUROC             #
# =============================================== #

  # Packages ----
  library("tidyverse")
  library("matrixStats")
  library("ggplot2")

  # General ----
  rm(list = ls())

  # Leave-one-study-out prediction ----
  df1 <- NULL
  for(s in 1:5){
    auc_tmp <- NULL
    for(a in c(1:100)){
      load(paste0("./CRC/AUROC/auroc.", as.character(s), ".", as.character(a), ".Rdata"))
      auc_tmp <- rbind(auc_tmp, AUROC.all)
    }
    df1 <- rbind(df1, auc_tmp %>% group_by(taxa.num, method, type) %>%
      summarize(AUC = mean(AUC, na.rm=TRUE)) %>%
        mutate(study = paste0("CRC",s), Method = method, linetype = type))
  }
  df1$taxa.num <- factor(df1$taxa.num, levels = unique(df1$taxa.num))
  df1$Method <- factor(df1$Method, levels = c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"))

  # Prepare plotting ----
  p1 <-  df1 %>% filter(study == "CRC1") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red", "purple" ,"orange", "grey", "green", "blue"))  +
    scale_linetype_manual(
      breaks = c("Original", "Study-effect corrected"),
      values = c("solid", "dashed"))  + theme_bw() +
    scale_y_continuous(limits = c(0.76, 0.84), breaks = seq(0.76, 0.84, by = 0.02)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size=14),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.box="vertical",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm'),
          strip.text = element_text(size = 18)) +
    facet_wrap(facets = vars(study)) + guides(
      color = guide_legend(order = 1, nrow = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Data")

  p2 <-  df1 %>% filter(study == "CRC2") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red", "purple" ,"orange", "grey", "green", "blue"))  +
    scale_linetype_manual(
      breaks = c("Original", "Study-effect corrected"),
      values = c("solid", "dashed"))  + theme_bw() +
    scale_y_continuous(limits = c(0.82, 0.9), breaks = seq(0.82, 0.9, by = 0.02)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size=14),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.box="vertical",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm'),
          strip.text = element_text(size = 18)) +
    facet_wrap(facets = vars(study)) + guides(
      color = guide_legend(order = 1, nrow = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Data")

  p3 <-  df1 %>% filter(study == "CRC3") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red", "purple" ,"orange", "grey", "green", "blue"))  +
    scale_linetype_manual(
      breaks = c("Original", "Study-effect corrected"),
      values = c("solid", "dashed"))  + theme_bw() +
    scale_y_continuous(limits = c(0.76, 0.84), breaks = seq(0.76, 0.84, by = 0.02)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size=14),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.box="vertical",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm'),
          strip.text = element_text(size = 18)) +
    facet_wrap(facets = vars(study)) + guides(
      color = guide_legend(order = 1, nrow = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Data")

  p4 <- df1 %>% filter(study == "CRC4") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red", "purple" ,"orange", "grey", "green", "blue"))  +
    scale_linetype_manual(
      breaks = c("Original", "Study-effect corrected"),
      values = c("solid", "dashed"))  + theme_bw() +
    scale_y_continuous(limits = c(0.73,0.812), breaks = seq(0.73,0.812, by = 0.02)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size=14),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.box="vertical",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm'),
          strip.text = element_text(size = 18)) +
    facet_wrap(facets = vars(study)) + guides(
      color = guide_legend(order = 1, nrow = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Data")

  p5 <- df1 %>% filter(study == "CRC5") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red", "purple" ,"orange", "grey", "green", "blue"))  +
    scale_linetype_manual(
      breaks = c("Original", "Study-effect corrected"),
      values = c("solid", "dashed"))  + theme_bw() +
    scale_y_continuous(limits = c(0.67,0.75), breaks = seq(0.67,0.75, by = 0.02)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size=14),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.box="vertical",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm'),
          strip.text = element_text(size = 18)) +
    facet_wrap(facets = vars(study)) + guides(
      color = guide_legend(order = 1, nrow = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Data")

  p6 <- df1 %>% group_by(taxa.num, Method, linetype) %>%  summarise(AUC = mean(AUC)) %>% add_column(study = "Combined") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "MMUPHin","CLR-LASSO","BW", "ALDEx2", "ANCOM-BC2"),
      values = c("red", "purple" ,"orange", "grey", "green", "blue"))  +
    scale_linetype_manual(
      breaks = c("Original", "Study-effect corrected"),
      values = c("solid", "dashed"))  + theme_bw() +
    scale_y_continuous(limits = c(0.74,0.82), breaks = seq(0.74,0.82, by = 0.02)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size=14),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.box="vertical",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm'),
          strip.text = element_text(size = 18)) +
    facet_wrap(facets = vars(study)) + guides(
      color = guide_legend(order = 1, nrow = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Data")

  # Generate figures ----
  ggp1 <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3,
                            common.legend = TRUE, legend = "bottom")

  pdf("./figures/Fig5.pdf", width = 14, height = 8, bg = "white")

  ggpubr::annotate_figure(ggp1, left = grid::textGrob("AUROC", rot = 90, vjust = 0.5,
                                                      hjust = -0.3, gp = grid::gpar(cex = 1.5)))

  dev.off()

