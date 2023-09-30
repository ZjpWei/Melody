  ##############################################
#                                                #
#    Figure 2 (e1) Stability: All v.s. LOSO      #
#                                                #
  ##############################################

  library("ggplot2")
  library("tidyverse")
  library("glmnet")

  rm(list = ls())
  # General
  data.all.loc <- "./CRC_meta_analysis/CRC_all_K401/"
  data.loso.loc <- "./CRC_meta_analysis/CRC_loso_K401/"
  load(paste0(data.all.loc, "prepare_data/data.rel.all.Rdata"))
  L <- length(data.rel)
  taxa.names <- colnames(data.rel[[1]]$Y)
  fdr <- 0.01

  # Original data
  ############################### Melody ####################################
  Jaccard <- c()
  methods <- c()
  study <- c()
  load(paste0(data.all.loc, "Models_original/Melody.model.all.Rdata"))
  global.meta <- list(names(which(Melody.model$coef!=0)))
  meta.rk <- list()
  for(l in 1:5){
    load(paste0(data.loso.loc, "Models_original/Melody.model.", as.character(l),".Rdata"))
    meta.rk[[l]] <- names(which(Melody.model$coef!=0))
  }
  meta.pr <- list()
  for(l in 1:5){
    meta.pr[[l]] <- length(intersect(meta.rk[[l]], global.meta[[1]]))/
      length(unique(c(meta.rk[[l]], global.meta[[1]])))
  }
  for(l in 1:L){
    Jaccard <- c(Jaccard, meta.pr[[l]])
    methods <- c(methods, "Melody")
    study <- c(study, as.character(l))
  }
  ############################### ANCOM-BC ##################################
  load(paste0(data.all.loc, "Models_original/ANCOMBC.model.all.Rdata"))
  global <- list(which(ANCOMBC.model$res$q_val$labels1 <= fdr))
  ANCOMBC.rk <- list()
  for(l in 1:5){
    load(paste0(data.loso.loc, "Models_original/ANCOMBC.model.", as.character(l),".Rdata"))
    ANCOMBC.rk[[l]] <- which(ANCOMBC.model$res$q_val$labels1 <= fdr)
  }
  ANCOMBC.pr <- list()
  for(l in 1:5){
    ANCOMBC.pr[[l]] <- length(intersect(ANCOMBC.rk[[l]], global[[1]]))/
      length(unique(c(ANCOMBC.rk[[l]], global[[1]])))
  }
  for(l in 1:L){
    Jaccard <- c(Jaccard, ANCOMBC.pr[[l]])
    methods <- c(methods, "ANCOM-BC")
    study <- c(study, as.character(l))
  }
  ############################### ALDEx2 ####################################
  load(paste0(data.all.loc, "Models_original/Aldex2.model.all.Rdata"))
  Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
  global <- list(which(Aldex2.fdr <= fdr))
  Aldex2.rk <- list()
  for(l in 1:5){
    load(paste0(data.loso.loc, "Models_original/Aldex2.model.", as.character(l), ".Rdata"))
    Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
    Aldex2.rk[[l]] <- which(Aldex2.fdr <= fdr)
  }
  Aldex2.pr <- list()
  for(l in 1:5){
    Aldex2.pr[[l]] <- length(intersect(Aldex2.rk[[l]], global[[1]]))/
      length(unique(c(Aldex2.rk[[l]], global[[1]])))
  }
  for(l in 1:L){
    Jaccard <- c(Jaccard, Aldex2.pr[[l]])
    methods <- c(methods, "ALDEx2")
    study <- c(study, as.character(l))
  }
  ############################### BW-prop ###################################
  load(paste0(data.all.loc, "Models_original/BW.prop.model.all.Rdata"))
  global <- list(which(BW.prop.model$q.val <= fdr))
  BW.rk <- list()
  for(l in 1:5){
    load(paste0(data.loso.loc, "Models_original/BW.prop.model.", as.character(l),".Rdata"))
    BW.rk[[l]] <- which(BW.prop.model$q.val <= fdr)
  }
  BW.pr <- list()
  for(l in 1:5){
    BW.pr[[l]] <- length(intersect(BW.rk[[l]], global[[1]]))/length(unique(c(BW.rk[[l]], global[[1]])))
  }
  for(l in 1:L){
    Jaccard <- c(Jaccard, BW.pr[[l]])
    methods <- c(methods, "BW")
    study <- c(study, as.character(l))
  }
  ############################### CLR-LASSO #################################
  load(paste0(data.all.loc, "Models_original/CLR.lasso.model.all.Rdata"))
  global.clrlasso <- list(names(which(clrlasso$glmnet.fit$beta[,clrlasso$lambda == clrlasso$lambda.min]!=0)))
  clrlasso.rk <- list()
  for(l in 1:5){
    load(paste0(data.loso.loc, "Models_original/CLR.lasso.model.", as.character(l), ".Rdata"))
    clrlasso.rk[[l]] <- names(which(clrlasso$glmnet.fit$beta[,clrlasso$lambda == clrlasso$lambda.min]!=0))
  }
  clrlasso.pr <- list()
  for(l in 1:5){
    clrlasso.pr[[l]] <- length(intersect(clrlasso.rk[[l]], global.clrlasso[[1]]))/
      length(unique(c(clrlasso.rk[[l]], global.clrlasso[[1]])))
  }

  for(l in 1:L){
    Jaccard <- c(Jaccard, clrlasso.pr[[l]])
    methods <- c(methods, "CLR-LASSO")
    study <- c(study, as.character(l))
  }
  df1 <- data.frame(Jaccard = Jaccard, methods = methods, study = study, type = "Original")

  # Batch-corrected data
  ############################### ANCOM-BC ##################################
  Jaccard <- c()
  methods <- c()
  study <- c()
  load(paste0(data.all.loc, "Models_batch_corrected/ANCOMBC.model.all.Rdata"))
  global <- list(which(ANCOMBC.model$res$q_val$labels1 <= fdr))
  ANCOMBC.rk <- list()
  for(l in 1:5){
    load(paste0(data.loso.loc, "Models_batch_corrected/ANCOMBC.model.", as.character(l),".Rdata"))
    ANCOMBC.rk[[l]] <- which(ANCOMBC.model$res$q_val$labels1 <= fdr)
  }
  ANCOMBC.pr <- list()
  for(l in 1:5){
    ANCOMBC.pr[[l]] <- length(intersect(ANCOMBC.rk[[l]], global[[1]]))/length(unique(c(ANCOMBC.rk[[l]], global[[1]])))
  }
  for(l in 1:L){
    Jaccard <- c(Jaccard, ANCOMBC.pr[[l]])
    methods <- c(methods, "ANCOM-BC")
    study <- c(study, as.character(l))
  }
  ############################### ALDEx2 ###################################
  load(paste0(data.all.loc, "Models_batch_corrected/Aldex2.model.all.Rdata"))
  Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
  global <- list(which(Aldex2.fdr <= fdr))
  Aldex2.rk <- list()
  for(l in 1:5){
    load(paste0(data.loso.loc, "Models_batch_corrected/Aldex2.model.", as.character(l), ".Rdata"))
    Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
    Aldex2.rk[[l]] <- which(Aldex2.fdr <= fdr)
  }
  Aldex2.pr <- list()
  for(l in 1:5){
    Aldex2.pr[[l]] <- length(intersect(Aldex2.rk[[l]], global[[1]]))/
      length(unique(c(Aldex2.rk[[l]], global[[1]])))
  }
  for(l in 1:L){
    Jaccard <- c(Jaccard, Aldex2.pr[[l]])
    methods <- c(methods, "ALDEx2")
    study <- c(study, as.character(l))
  }
  ################################ BW-prop ##################################
  load(paste0(data.all.loc, "Models_batch_corrected/BW.prop.model.all.Rdata"))
  global <- list(which(BW.prop.model$q.val <= fdr))
  BW.rk <- list()
  for(l in 1:5){
    load(paste0(data.loso.loc, "Models_batch_corrected/BW.prop.model.", as.character(l),".Rdata"))
    BW.rk[[l]] <- which(BW.prop.model$q.val <= fdr)
  }
  BW.pr <- list()
  for(l in 1:5){
    BW.pr[[l]] <- length(intersect(BW.rk[[l]], global[[1]]))/
      length(unique(c(BW.rk[[l]], global[[1]])))
  }
  for(l in 1:L){
    Jaccard <- c(Jaccard, BW.pr[[l]])
    methods <- c(methods, "BW")
    study <- c(study, as.character(l))
  }
  ############################### CLR-LASSO #################################
  load(paste0(data.all.loc, "Models_batch_corrected/CLR.lasso.model.all.Rdata"))
  global.clrlasso <- list(names(which(clrlasso$glmnet.fit$beta[,clrlasso$lambda == clrlasso$lambda.min]!=0)))
  clrlasso.rk <- list()
  for(l in 1:5){
    load(paste0(data.loso.loc, "Models_batch_corrected/CLR.lasso.model.", as.character(l), ".Rdata"))
    clrlasso.rk[[l]] <- names(which(clrlasso$glmnet.fit$beta[,clrlasso$lambda == clrlasso$lambda.min]!=0))
  }
  clrlasso.pr <- list()
  for(l in 1:5){
    clrlasso.pr[[l]] <- length(intersect(clrlasso.rk[[l]], global.clrlasso[[1]]))/
      length(unique(c(clrlasso.rk[[l]], global.clrlasso[[1]])))
  }
  for(l in 1:L){
    Jaccard <- c(Jaccard, clrlasso.pr[[l]])
    methods <- c(methods, "CLR-LASSO")
    study <- c(study, as.character(l))
  }
  df2 <- data.frame(Jaccard = Jaccard, methods = methods, study = study, type = "Batch-corrected")
  df <- rbind(df1, df2)
  df$methods <- factor(df$methods, c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"))
  df$type <- factor(df$type, c("Original", "Batch-corrected"))

  # Prepare plotting
  Jac_maxs <- c()
  Jac_mins <- c()
  for(l in 1:nrow(df)){
    Jac_maxs <- c(Jac_maxs, max(df %>% filter(methods == df$methods[l]) %>%
                                  filter(type == df$type[l]) %>%
                                  pull(Jaccard)))
    Jac_mins <- c(Jac_mins, min(df %>% filter(methods == df$methods[l]) %>%
                                  filter(type == df$type[l]) %>%
                                  pull(Jaccard)))
  }
  df$Jac_max <- Jac_maxs
  df$Jac_min <- Jac_mins

  # Customizing the output
  pdf("./figures/figure2_e1.pdf",         # File name
      width = 4.72, height = 4.45, # Width and height in inches
      bg = "white")

  # Creating a plot
  df  %>% ggplot(aes(x=methods, y=Jaccard, group = type, color = methods)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) + ylim(0,1) +
    geom_errorbar(
      aes(x = methods,ymin = Jac_min, ymax = Jac_max, color= methods,linetype = type),
      width = 0, position = position_dodge(0.5)) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" )
    ) +  theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.title = element_blank(),
               axis.text.y = element_text(size = 15),
               axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
               panel.background = element_rect(fill = 'white'),
               panel.border = element_rect(colour = "black", fill=NA),
               legend.position = "none",  text = element_text(size=15))

  # Closing the graphical device
  dev.off()
