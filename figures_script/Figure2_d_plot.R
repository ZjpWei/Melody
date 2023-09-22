  ##############################################
#                                                #
#         Figure 2 (d) random split and          #
#         Leave-on-study-out Prediction          #
#                                                #
  ##############################################

  library("tidyverse")
  library("matrixStats")
  library("randomForest")
  library("ggplot2")
  
  rm(list = ls())
  # Random split prediction
  data.loc <- "./Prediction/CRC_ltcvo_K401/"
  
  # Prepare plotting
  taxa.num <- NULL
  method <- NULL
  linetype <- NULL
  AUC <- NULL
  for(methodology in c("Melody", "clrlasso", "Aldex2", "ANCOM", "BW.p")){
    ## Defining the AUROC matrix for original data
    auc_tmp <- NULL 
    for(a in c(1:100))
    {
      if(file.exists(paste0(data.loc, "AUROC/auroc.", as.character(a), ".prop.",methodology , ".Rdata"))){
        load(paste0(data.loc, "AUROC/auroc.", as.character(a), ".prop.",methodology , ".Rdata"))
        auc_tmp <- rbind(auc_tmp, auroc.all$AUC)
      }
    }
    AUC <- c(AUC, colMeans(auc_tmp, na.rm = TRUE))
    taxa.num <- c(taxa.num, auroc.all$taxa.num)
    method <- c(method, auroc.all$method)
    linetype <- c(linetype, rep("Original", length(auroc.all$method)))
    # Defining the AUROC matrix for Batch-corrected data
    if(methodology %in% c("clrlasso", "Aldex2", "ANCOM", "BW.p")){
      auc_tmp <- NULL
      for(a in c(1:100))
      {
        if(file.exists(paste0(data.loc, "AUROC_batch/auroc.", as.character(a), ".prop.",methodology , ".Rdata"))){
          load(paste0(data.loc, "AUROC_batch/auroc.", as.character(a), ".prop.",methodology , ".Rdata"))
          auc_tmp <- rbind(auc_tmp, auroc.all$AUC)
        }
      }
      AUC <- c(AUC, colMeans(auc_tmp, na.rm = TRUE))
      taxa.num <- c(taxa.num, auroc.all$taxa.num)
      method <- c(method, auroc.all$method)
      linetype <- c(linetype, rep("Batch-corrected", length(auroc.all$method)))
    }
  }
  df1 <- data.frame(AUC = AUC, method = method, taxa.num = taxa.num, linetype = linetype)
  df1$taxa.num <- factor(df1$taxa.num, levels = unique(df1$taxa.num))
  df1$method[df1$method == "Melody_1"] <- "Melody"
  df1$method[df1$method == "Aldex2"] <- "ALDEx2"
  df1$method[df1$method == "ANCOM"] <- "ANCOM-BC"
  df1$method[df1$method == "clrlasso"] <- "CLR-LASSO"
  df1$method[df1$method == "BW.p"] <- "BW"
  df1$method <- factor(df1$method, levels = c("Melody","ALDEx2","ANCOM-BC","BW","CLR-LASSO"))
  
  # Customizing the output
  pdf("./figures/figure2_d_RS.pdf",         # File name
      width = 7.07, height = 5.25, # Width and height in inches
      bg = "white")
  
  ggplot(df1, aes(x=taxa.num, y=AUC, group = interaction(method,linetype), linetype = linetype)) + 
    ylim(0.75,0.85) +
    geom_line(aes(color = method)) + 
    geom_point(aes(color = method), pch = 18) +
    ggtitle("Random split") +
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
          legend.key.height = unit(0.7, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm')) + 
    labs(linetype="Pooled data")
  
  # Closing the graphical device
  dev.off()
 
  ##############################################################################
  rm(list = ls())
  # Leave-one-study-out prediction
  data.loc <- "./Prediction/CRC_locvo_K401/"
  taxa.num <- NULL
  method <- NULL
  linetype <- NULL
  study <- NULL
  AUC <- NULL
  for(methodology in c("Melody_1", "clrlasso","Aldex2", "ANCOM", "BW.p")){
    for(s in 1:5){
      ## Defining the AUROC matrix for ANCOM BC
      auc_tmp <- NULL 
      for(a in c(1:100))
      {
        load(paste0(data.loc, "AUROC/auroc.", as.character(s), ".", as.character(a), ".prop.",methodology , ".Rdata"))
        auc_tmp <- rbind(auc_tmp, auroc.all$AUC)
      }
      AUC <- c(AUC, colMeans(auc_tmp, na.rm = TRUE))
      study <- c(study, auroc.all$s)
      taxa.num <- c(taxa.num, auroc.all$taxa.num)
      method <- c(method, auroc.all$method)
      linetype <- c(linetype, rep("Original", length(auroc.all$method)))
    }
    if(methodology != "Melody_1"){
      for(s in 1:5){
        ## Defining the AUROC matrix for ANCOM BC
        auc_tmp <- NULL 
        for(a in c(1:100)){
          loc <- paste0(data.loc, "AUROC_batch/auroc.", as.character(s), ".", as.character(a), ".prop.",methodology , ".Rdata")
          if(file.exists(loc)){
            load(loc)
            auc_tmp <- rbind(auc_tmp, auroc.all$AUC)
          }
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
  df1$method[df1$method == "Melody_1"] <- "Melody"
  df1$method[df1$method == "Aldex2"] <- "ALDEx2"
  df1$method[df1$method == "ANCOM"] <- "ANCOM-BC"
  df1$method[df1$method == "clrlasso"] <- "CLR-LASSO"
  df1$method[df1$method == "BW.p"] <- "BW"
  df1$method <- factor(df1$method, levels = c("Melody","ALDEx2","ANCOM-BC","BW","CLR-LASSO"))
  
  pdf("./figures/figure2_d_LOSO.pdf",         # File name
      width = 7.07, height = 5.25, # Width and height in inches
      bg = "white")
  
  df1 %>% group_by(taxa.num, method, linetype) %>%  summarise(AUC = mean(AUC)) %>%
    ggplot(aes(x=taxa.num, y=AUC, group = interaction(method,linetype), linetype = linetype)) +
    geom_line(aes(color = method)) + ylim(0.725,0.825) +
    geom_point(aes(color = method), pch = 18) +
    ggtitle("LOSO") + 
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
          legend.key.width = unit(1, 'cm')) + 
    labs(linetype="Pooled data", color="")
  
  dev.off()