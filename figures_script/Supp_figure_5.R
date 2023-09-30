  ##############################################
#                                                #
#             Supplementary figure 5             #
#                                                #
  ##############################################

  library("tidyverse")
  library("matrixStats")
  library("ggplot2")

  # Prepare plotting
  rm(list = ls())
  # Leave-one-study-out prediction
  data.loc <- "./CRC_meta_analysis/Prediction/LOSO/"
  taxa.num <- NULL
  method <- NULL
  linetype <- NULL
  study <- NULL
  AUC <- NULL
  for(methodology in c("Melody", "clrlasso","Aldex2", "ANCOM", "BW")){
    for(s in 1:5){
      ## Defining the AUROC matrix for ANCOM BC
      auc_tmp <- NULL
      for(a in c(1:100)){
        load(paste0(data.loc, "AUORC_original/auroc.", as.character(s), ".", as.character(a), ".prop.",methodology , ".Rdata"))
        auc_tmp <- rbind(auc_tmp, auroc.all$AUC)
      }
      AUC <- c(AUC, colMeans(auc_tmp, na.rm = TRUE))
      study <- c(study, auroc.all$s)
      taxa.num <- c(taxa.num, auroc.all$taxa.num)
      method <- c(method, auroc.all$method)
      linetype <- c(linetype, rep("Original", length(auroc.all$method)))
    }
    if(methodology %in% c("clrlasso", "Aldex2", "ANCOM", "BW")){
      for(s in 1:5){
        ## Defining the AUROC matrix for ANCOM BC
        auc_tmp <- NULL
        for(a in c(1:100)){
          load(paste0(data.loc, "AUROC_batch_corrected/auroc.", as.character(s), ".", as.character(a), ".prop.",methodology , ".Rdata"))
          auc_tmp <- rbind(auc_tmp, auroc.all$AUC)
        }
        AUC <- c(AUC, colMeans(auc_tmp, na.rm = TRUE))
        study <- c(study, auroc.all$s)
        taxa.num <- c(taxa.num, as.character((1:16) * 5))
        method <- c(method, auroc.all$method)
        linetype <- c(linetype, rep("Batch-corrected", length(auroc.all$method)))
      }
    }
  }
  df1 <- data.frame(AUC = AUC, method = method, taxa.num = taxa.num, linetype = linetype, study = study)
  df1$taxa.num <- factor(df1$taxa.num, levels = unique(df1$taxa.num))
  df1$method[df1$method == "Aldex2"] <- "ALDEx2"
  df1$method[df1$method == "ANCOM"] <- "ANCOM-BC"
  df1$method[df1$method == "clrlasso"] <- "CLR-LASSO"
  df1$method <- factor(df1$method, levels = c("Melody","ALDEx2","ANCOM-BC","BW","CLR-LASSO"))

  # Plotting
  p1 <-  df1 %>% filter(study == "1") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(method,linetype), linetype = linetype)) +
    geom_line(aes(color = method)) + ylim(0.7,0.85) +
    geom_point(aes(color = method), pch = 18) +
    ggtitle("CRC1") +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
      values = c("solid", "dashed"))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size=17),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'))

  p2 <-  df1 %>% filter(study == "2") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(method,linetype), linetype = linetype)) +
    geom_line(aes(color = method)) + ylim(0.75, 0.905) +
    geom_point(aes(color = method), pch = 18) +
    ggtitle("CRC2") +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
      values = c("solid", "dashed"))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size=17),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'))

  p3 <-  df1 %>% filter(study == "3") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(method,linetype), linetype = linetype)) +
    geom_line(aes(color = method)) + ylim(0.75, 0.85) +
    geom_point(aes(color = method), pch = 18) +
    ggtitle("CRC3") +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
      values = c("solid", "dashed"))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size=17),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'))

  ## b
  p4 <- df1 %>% filter(study == "4") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(method,linetype), linetype = linetype)) +
    geom_line(aes(color = method)) + ylim(0.70,0.82) +
    geom_point(aes(color = method), pch = 18) +
    ggtitle("CRC4") +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
      values = c("solid", "dashed"))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size=17),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'))
  p5 <-  df1 %>% filter(study == "5") %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(method,linetype), linetype = linetype)) +
    geom_line(aes(color = method)) + ylim(0.65,0.75) +
    geom_point(aes(color = method), pch = 18) +
    ggtitle("CRC5") +
    xlab("Top"~italic("N")~"selected speices") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" ))  +
    scale_linetype_manual(
      breaks = c("Original", "Batch-corrected"),
      values = c("solid", "dashed"))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size=17),
          plot.title = element_text(size=20,hjust = 0.5, vjust = -3),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'))
  ggp1 <- ggpubr::ggarrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3)

  pdf("./figures/Supp_5.pdf",         # File name
      width = 20, height = 10, # Width and height in inches
      bg = "white")

  ggp1


  dev.off()

