# =============================================== #
#              Figure 5: Stability                #
# =============================================== #
  
  # Packages ----
  library("ggplot2")
  library("tidyverse")
  library("glmnet")
  library("ggtext")
  
  # General ----
  rm(list = ls())
  fdr.cut <- 0.05
  method.name <- c("Melody", "ALDEx2", "ANCOMBC2", "MMUPHin", "clrlasso", "BW")
  
  # K401 v.s. K849 ----
  Jac.dist <- NULL
  methods <- NULL
  data.type <- NULL
  for(ss in c(1:100)){
    load(paste0("./CRC/Stability/Model.summary.K401.all.",ss,".Rdata"))
    assign("res_K401", Model.summary)
    load(paste0("./CRC/Stability/Model.summary.K849.all.",ss,".Rdata"))
    assign("res_K849", Model.summary)
    for(meth in method.name){
      if(meth == "Melody"){
        res_tmp <- res_K401 %>% dplyr::transmute(feature, coef_K401 = Melody_coef) %>%
          dplyr::left_join(res_K849 %>% dplyr::transmute(feature, coef_K849 = Melody_coef),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(coef_K401 != 0) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(coef_K849 != 0) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Original")
      }else if(meth == "clrlasso"){
        ## Original data
        res_tmp <- res_K401 %>% dplyr::transmute(feature, coef_K401 = clrlasso_org_coef) %>%
          dplyr::left_join(res_K849 %>% dplyr::transmute(feature, coef_K849 = clrlasso_org_coef),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(coef_K401 != 0) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(coef_K849 != 0) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Original")
        
        ## Batch-corrected data
        res_tmp <- res_K401 %>% dplyr::transmute(feature, coef_K401 = clrlasso_bat_coef) %>%
          dplyr::left_join(res_K849 %>% dplyr::transmute(feature, coef_K849 = clrlasso_bat_coef),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(coef_K401 != 0) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(coef_K849 != 0) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Study-effect corrected")
      }else{
        ## Original data
        res_tmp <- res_K401 %>% dplyr::transmute(feature, fdr_K401 = get(paste0(meth, "_org_fdr"))) %>%
          dplyr::left_join(res_K849 %>% dplyr::transmute(feature, fdr_K849 = get(paste0(meth, "_org_fdr"))),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(fdr_K401 <= fdr.cut) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(fdr_K849 <= fdr.cut) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Original")
        
        ## Batch-corrected data
        res_tmp <- res_K401 %>% dplyr::transmute(feature, fdr_K401 = get(paste0(meth, "_bat_fdr"))) %>%
          dplyr::left_join(res_K849 %>% dplyr::transmute(feature, fdr_K849 = get(paste0(meth, "_bat_fdr"))),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(fdr_K401 <= fdr.cut) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(fdr_K849 <= fdr.cut) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Study-effect corrected")
      }
    }
  }
  df.K401.K849 <- data.frame(Jaccard = Jac.dist, methods = methods, type = data.type) %>% 
    dplyr::transmute(methods, Jaccard, type = type)
  df.K401.K849$methods[df.K401.K849$methods == "ANCOMBC2"] <- "ANCOM-BC2"
  df.K401.K849$methods[df.K401.K849$methods == "clrlasso"] <- "CLR-LASSO"
  df.K401.K849$methods <- factor(df.K401.K849$methods, levels =c("Melody","MMUPHin", "CLR-LASSO","BW","ALDEx2","ANCOM-BC2"))
  df.K401.K849$type <- factor(df.K401.K849$type, levels = c("Original", "Study-effect corrected"))

  p1 <- df.K401.K849 %>% ggplot(aes(x=methods, y=Jaccard, color = methods, linetype = type)) +
    geom_boxplot(position = position_dodge2(preserve = "single"), width=0.4) + 
    ylim(0,1) + ylab("Jaccard Index") +
    scale_color_manual(
      breaks = c("Melody","MMUPHin", "CLR-LASSO","BW","ALDEx2", "ANCOM-BC2"),
      values = c("red", "purple", "orange","grey","green" ,"blue")) +
    ggtitle("<i>K</i> = 401 vs <i>K</i> = 849 species") + xlab("Study") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_markdown(hjust = 0.5, face = "bold", size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",  text = element_text(size=15)) +
    labs(linetype="Data") + guides(colour = guide_legend(nrow = 2))
  
  # K401 v.s. order ----
  Jac.dist <- NULL
  methods <- NULL
  data.type <- NULL
  for(ss in 1:100){
    load(paste0("./CRC/Stability/Model.summary.K401.all.",ss,".Rdata"))
    assign("res_K401", Model.summary)
    load(paste0("./CRC/Stability/Model.summary.order.all.",ss,".Rdata"))
    assign("res_order", Model.summary)
    for(meth in method.name){
      if(meth == "Melody"){
        res_tmp <- res_order %>% dplyr::transmute(feature, coef_order = Melody_coef) %>%
          dplyr::left_join(res_K401 %>% dplyr::transmute(feature, coef_K401 = Melody_coef),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(coef_K401 != 0) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(coef_order != 0) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Original")
      }else if(meth == "clrlasso"){
        ## Original data
        res_tmp <- res_order %>% dplyr::transmute(feature, coef_order = clrlasso_org_coef) %>%
          dplyr::left_join(res_K401 %>% dplyr::transmute(feature, coef_K401 = clrlasso_org_coef),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(coef_K401 != 0) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(coef_order != 0) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Original")
        
        ## Batch-corrected data
        res_tmp <- res_order %>% dplyr::transmute(feature, coef_order = clrlasso_bat_coef) %>%
          dplyr::left_join(res_K401 %>% dplyr::transmute(feature, coef_K401 = clrlasso_bat_coef),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(coef_K401 != 0) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(coef_order != 0) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Study-effect corrected")
      }else{
        ## Original data
        res_tmp <- res_order %>% dplyr::transmute(feature, fdr_order = get(paste0(meth, "_org_fdr"))) %>%
          dplyr::left_join(res_K401 %>% dplyr::transmute(feature, fdr_K401 = get(paste0(meth, "_org_fdr"))),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(fdr_K401 <= fdr.cut) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(fdr_order <= fdr.cut) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Original")
        
        ## Batch-corrected data
        res_tmp <- res_order %>% dplyr::transmute(feature, fdr_order = get(paste0(meth, "_bat_fdr"))) %>%
          dplyr::left_join(res_K401 %>% dplyr::transmute(feature, fdr_K401 = get(paste0(meth, "_bat_fdr"))),
                           by = "feature")
        feature.1 <- res_tmp %>% dplyr::filter(fdr_K401 <= fdr.cut) %>% dplyr::pull(feature)
        feature.2 <- res_tmp %>% dplyr::filter(fdr_order <= fdr.cut) %>% dplyr::pull(feature)
        Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
        methods <- c(methods, meth)
        data.type <- c(data.type, "Study-effect corrected")
      }
    }
  }
  df.K401.order <- data.frame(Jaccard = Jac.dist, methods = methods, type = data.type) %>% 
    dplyr::transmute(methods, Jaccard, type = type)
  df.K401.order$methods[df.K401.order$methods == "ANCOMBC2"] <- "ANCOM-BC2"
  df.K401.order$methods[df.K401.order$methods == "clrlasso"] <- "CLR-LASSO"
  df.K401.order$methods <- factor(df.K401.order$methods, levels = c("Melody","MMUPHin", "CLR-LASSO","BW","ALDEx2","ANCOM-BC2"))
  df.K401.order$type <- factor(df.K401.order$type, levels = c("Original", "Study-effect corrected"))
  
  p2 <- df.K401.order %>% ggplot(aes(x=methods, y=Jaccard, color = methods, linetype = type)) +
    geom_boxplot(position = position_dodge2(preserve = "single"), width=0.4) + 
    ylim(0,1) + ylab("Jaccard Index") +
    scale_color_manual(
      breaks = c("Melody","MMUPHin", "CLR-LASSO","BW","ALDEx2","ANCOM-BC2"),
      values = c("red", "purple", "orange","grey","green" ,"blue")) +
    ggtitle("<i>K</i> = 401 vs <i>K</i> = 267 species") + xlab("Study") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_markdown(hjust = 0.5, face = "bold", size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",  text = element_text(size=15)) +
    labs(linetype="Data") + guides(colour = guide_legend(nrow = 2))
  
  # All v.s. loso ----
  df.K401.all <- NULL
  for(s in 1:5){
    Jac.dist <- NULL
    methods <- NULL
    data.type <- NULL
    for(ss in c(1:100)){
      load(paste0("./CRC/Stability/Model.summary.K401.all.",ss,".Rdata"))
      assign("res_all", Model.summary)
      load(paste0("./CRC/Stability/Model.summary.K401.",s,".",ss,".Rdata"))
      assign("res_loso", Model.summary)
      for(meth in method.name){
        if(meth == "Melody"){
          res_tmp <- res_all %>% dplyr::transmute(feature, coef_all = Melody_coef) %>%
            dplyr::left_join(res_loso %>% dplyr::transmute(feature, coef_loso = Melody_coef),
                             by = "feature")
          feature.1 <- res_tmp %>% dplyr::filter(coef_all != 0) %>% dplyr::pull(feature)
          feature.2 <- res_tmp %>% dplyr::filter(coef_loso != 0) %>% dplyr::pull(feature)
          Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
          methods <- c(methods, meth)
          data.type <- c(data.type, "Original")
        }else if(meth == "clrlasso"){
          ## Original data
          res_tmp <- res_all %>% dplyr::transmute(feature, coef_all = clrlasso_org_coef) %>%
            dplyr::left_join(res_loso %>% dplyr::transmute(feature, coef_loso = clrlasso_org_coef),
                             by = "feature")
          feature.1 <- res_tmp %>% dplyr::filter(coef_all != 0) %>% dplyr::pull(feature)
          feature.2 <- res_tmp %>% dplyr::filter(coef_loso != 0) %>% dplyr::pull(feature)
          Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
          methods <- c(methods, meth)
          data.type <- c(data.type, "Original")
          
          ## Batch-corrected data
          res_tmp <- res_all %>% dplyr::transmute(feature, coef_all = clrlasso_bat_coef) %>%
            dplyr::left_join(res_loso %>% dplyr::transmute(feature, coef_loso = clrlasso_bat_coef),
                             by = "feature")
          feature.1 <- res_tmp %>% dplyr::filter(coef_all != 0) %>% dplyr::pull(feature)
          feature.2 <- res_tmp %>% dplyr::filter(coef_loso != 0) %>% dplyr::pull(feature)
          Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
          methods <- c(methods, meth)
          data.type <- c(data.type, "Study-effect corrected")
        }else{
          ## Original data
          res_tmp <- res_all %>% dplyr::transmute(feature, fdr_all = get(paste0(meth, "_org_fdr"))) %>%
            dplyr::left_join(res_loso %>% dplyr::transmute(feature, fdr_loso = get(paste0(meth, "_org_fdr"))),
                             by = "feature")
          feature.1 <- res_tmp %>% dplyr::filter(fdr_all <= fdr.cut) %>% dplyr::pull(feature)
          feature.2 <- res_tmp %>% dplyr::filter(fdr_loso <= fdr.cut) %>% dplyr::pull(feature)
          Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
          methods <- c(methods, meth)
          data.type <- c(data.type, "Original")
          
          ## Batch-corrected data
          res_tmp <- res_all %>% dplyr::transmute(feature, fdr_all = get(paste0(meth, "_bat_fdr"))) %>%
            dplyr::left_join(res_loso %>% dplyr::transmute(feature, fdr_loso = get(paste0(meth, "_bat_fdr"))),
                             by = "feature")
          feature.1 <- res_tmp %>% dplyr::filter(fdr_all <= fdr.cut) %>% dplyr::pull(feature)
          feature.2 <- res_tmp %>% dplyr::filter(fdr_loso <= fdr.cut) %>% dplyr::pull(feature)
          Jac.dist <- c(Jac.dist, length(intersect(feature.1, feature.2)) / length(unique(c(feature.1, feature.2))))
          methods <- c(methods, meth)
          data.type <- c(data.type, "Study-effect corrected")
        }
      }
    }
    df.tmp <- data.frame(Jaccard = Jac.dist, methods = methods, type = data.type) %>% 
      dplyr::transmute(methods, Jaccard, study = s, type = type)
    df.K401.all <- rbind(df.K401.all, df.tmp)
  }
  df.K401.all$methods[df.K401.all$methods == "ANCOMBC2"] <- "ANCOM-BC2"
  df.K401.all$methods[df.K401.all$methods == "clrlasso"] <- "CLR-LASSO"
  df.K401.all$methods <- factor(df.K401.all$methods, levels = c("Melody","MMUPHin", "CLR-LASSO","BW","ALDEx2","ANCOM-BC2"))
  df.K401.all$type <- factor(df.K401.all$type, levels = c("Original", "Study-effect corrected"))

  p3 <- df.K401.all %>% ggplot(aes(x=factor(study), y=Jaccard, color = methods, linetype = type)) +
    geom_boxplot(position = position_dodge2(preserve = "single"), width=0.6) + 
    ylim(0,1) + ylab("Jaccard Index") +
    scale_color_manual(
      breaks = c("Melody","MMUPHin", "CLR-LASSO","BW","ALDEx2","ANCOM-BC2"),
      values = c("red", "purple", "orange","grey","green" ,"blue")) +
    ggtitle("All vs LOSO data") + labs(x = "Left-out study") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          panel.background = element_rect(fill = 'white'),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "bottom",  text = element_text(size=15),
          legend.box="vertical",
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.background = element_rect(fill = 'NA'),
          legend.key =  element_rect(colour = "transparent", fill = "transparent"),
          legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'),
          strip.text = element_text(size = 14)) +
    labs(linetype="Data") + guides(colour = guide_legend(nrow = 1))

  
  pdf("./figures/Fig4_K401_K849.pdf", width = 4, height = 3.00, bg = "white")

  p1

  dev.off()

  pdf("./figures/Fig4_K401_order.pdf", width = 4, height = 3.00, bg = "white")

  p2

  dev.off()

  pdf("./figures/Fig4_All_loso.pdf", width = 10.82, height = 4.5, bg = "white")

  p3

  dev.off()
