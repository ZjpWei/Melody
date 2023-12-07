# =============================================== #
#                    Figure 6                     #
# =============================================== #

  # Packages ----
  library("tidyverse")
  library("matrixStats")
  library("ggplot2")

  # General ----
  rm(list = ls())
  
  # Leave-one-study-out prediction ----
  data.loc <- "./CRC_meta_analysis/Prediction/LOSO/"
  taxa.num <- NULL
  Method <- NULL
  linetype <- NULL
  study <- NULL
  AUC <- NULL
  for(Methodology in c("Melody", "clrlasso","Aldex2", "ANCOM", "BW")){
    for(s in 1:5){
      auc_tmp <- NULL
      for(a in c(1:100)){
        load(paste0(data.loc, "AUROC_original/auroc.", as.character(s), ".", as.character(a), ".prop.",Methodology , ".Rdata"))
        auc_tmp <- rbind(auc_tmp, auroc.all$AUC)
      }
      AUC <- c(AUC, colMeans(auc_tmp, na.rm = TRUE))
      study <- c(study, paste0("CRC", auroc.all$s))
      taxa.num <- c(taxa.num, auroc.all$taxa.num)
      Method <- c(Method, auroc.all$method)
      linetype <- c(linetype, rep("Original", length(auroc.all$method)))
    }
    if(Methodology %in% c("clrlasso", "Aldex2", "ANCOM", "BW")){
      for(s in 1:5){
        auc_tmp <- NULL
        for(a in c(1:100)){
          load(paste0(data.loc, "AUROC_batch_corrected/auroc.", as.character(s), ".", as.character(a), ".prop.",Methodology , ".Rdata"))
          auc_tmp <- rbind(auc_tmp, auroc.all$AUC)
        }
        AUC <- c(AUC, colMeans(auc_tmp, na.rm = TRUE))
        study <- c(study, paste0("CRC", auroc.all$s))
        taxa.num <- c(taxa.num, as.character((2:16) * 5))
        Method <- c(Method, auroc.all$method)
        linetype <- c(linetype, rep("Batch-corrected", length(auroc.all$method)))
      }
    }
  }
  df1 <- data.frame(AUC = AUC, Method = Method, taxa.num = taxa.num, linetype = linetype, study = study)
  df1$taxa.num <- factor(df1$taxa.num, levels = unique(df1$taxa.num))
  df1$Method[df1$Method == "Aldex2"] <- "ALDEx2"
  df1$Method[df1$Method == "ANCOM"] <- "ANCOM-BC"
  df1$Method[df1$Method == "clrlasso"] <- "CLR-LASSO"
  df1$Method <- factor(df1$Method, levels = c("Melody","ALDEx2","ANCOM-BC","BW","CLR-LASSO"))

  # Prepare plotting ----
  p1 <-  df1 %>% filter(study == "CRC1") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) + 
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
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
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Pooled data")

  p2 <-  df1 %>% filter(study == "CRC2") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) + 
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
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
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Pooled data")

  p3 <-  df1 %>% filter(study == "CRC3") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
      values = c("solid", "dashed"))  + theme_bw() +
    scale_y_continuous(limits = c(0.77, 0.85), breaks = seq(0.77, 0.85, by = 0.02)) +
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
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Pooled data")

  p4 <- df1 %>% filter(study == "CRC4") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) + 
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
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
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Pooled data")
  
  p5 <- df1 %>% filter(study == "CRC5") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) + 
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
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
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Pooled data")
  
  p6 <- df1 %>% group_by(taxa.num, Method, linetype) %>%  summarise(AUC = mean(AUC)) %>% add_column(study = "Combined") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(Method,linetype), linetype = linetype)) +
    geom_line(aes(color = Method)) +
    geom_point(aes(color = Method), pch = 18) +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
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
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    ) + labs(linetype="Pooled data")
    
  # Generate figures ----
  ggp1 <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3, 
                            common.legend = TRUE, legend = "bottom")

  pdf("./figures/figure6.pdf", width = 14, height = 8, bg = "white")

  ggpubr::annotate_figure(ggp1, left = grid::textGrob("AUROC", rot = 90, vjust = 0.5, 
                                                      hjust = -0.3, gp = grid::gpar(cex = 1.5)))
  
  dev.off()

