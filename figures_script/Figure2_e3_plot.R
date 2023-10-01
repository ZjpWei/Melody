  ##############################################
#                                                #
# Figure 2 (e3) Stability: K401 v.s. K267(order) #
#                                                #
  ##############################################

  library("ggplot2")
  library("tidyverse")
  library("glmnet")

  rm(list = ls())
  # General
  load("./CRC_meta_analysis/data/tax.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  taxa.name1 <- colnames(data.rel[[1]]$Y)
  load("./CRC_meta_analysis/CRC_all_order/prepare_data/data.rel.all.Rdata")
  taxa.name2 <- colnames(data.rel[[1]]$Y)
  filters <- match(taxa.name2, taxa.name1)
  tax.species <- gsub("[][]", "",unlist(regmatches(colnames(data.rel[[1]]$Y),
                      gregexpr("\\[.*?\\]",colnames(data.rel[[1]]$Y)))))
  tax.sub <- tax[match(tax.species, tax$OTUID),]
  fdr <- 0.01
  replicates <- c(1:5,"all")
  data.loc.1 <- "./CRC_meta_analysis/CRC_loso_K401/"
  data.loc.2 <- "./CRC_meta_analysis/CRC_loso_order/"

  # Original data
  ############################### Melody ####################################
  Jaccard <- NULL
  Ham_dist <- NULL
  methods <- NULL
  study <- NULL
  meta.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_order/Models_original/Melody.model.all.Rdata")
      melody.id <- names(Melody.model$coef)
    }else{
      load(paste0(data.loc.2, "Models_original/Melody.model.", as.character(l),".Rdata"))
      melody.id <- names(Melody.model$coef)
    }
    tmp.id <- which(Melody.model$coef != 0)
    meta.rk.1[[k]] <- names(tmp.id)
    k <- k + 1
  }
  meta.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/Melody.model.all.Rdata")
      melody.id.2 <- names(Melody.model$coef)
      filters.id <- intersect(melody.id, melody.id.2)
    }else{
      load(paste0(data.loc.1, "Models_original/Melody.model.", as.character(l),".Rdata"))
      melody.id.2 <- names(Melody.model$coef)
      filters.id <- intersect(melody.id, melody.id.2)
    }
    meta.rk.2[[k]] <- names(which(Melody.model$coef[filters.id] != 0))
    k <- k + 1
  }
  meta.pr <- list()
  for(l in 1:length(replicates)){
    percent.tx <- length(intersect(meta.rk.1[[l]], meta.rk.2[[l]]))/
      length(unique(c(meta.rk.1[[l]], meta.rk.2[[l]])))
    meta.pr[[l]] <- percent.tx
  }
  for(l in 1:length(replicates)){
    Jaccard <- c(Jaccard, meta.pr[[l]])
    methods <- c(methods, "Melody")
  }
  study <- c(rep("loso", 5), "all")
  ############################### ANCOM-BC ##################################
  ANCOMBC.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_order/Models_original/ANCOMBC.model.all.Rdata")
      ANCOMBC.rk.1[[k]] <- which(ANCOMBC.model$res$q_val$labels1 <= fdr)
    }else{
      load(paste0(data.loc.2, "Models_original/ANCOMBC.model.", as.character(l),".Rdata"))
      ANCOMBC.rk.1[[k]] <- which(ANCOMBC.model$res$q_val$labels1 <= fdr)
    }
    k <- k + 1
  }
  ANCOMBC.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/ANCOMBC.model.all.Rdata")
      ANCOMBC.rk.2[[k]] <- which(ANCOMBC.model$res$q_val$labels1[filters] <= fdr)
    }else{
      load(paste0(data.loc.1, "Models_original/ANCOMBC.model.", as.character(l),".Rdata"))
      ANCOMBC.rk.2[[k]] <- which(ANCOMBC.model$res$q_val$labels1[filters] <= fdr)
    }
    k <- k + 1
  }
  ANCOMBC.pr <- list()
  for(l in 1:length(replicates)){
    percent.tx <- length(intersect(ANCOMBC.rk.1[[l]], ANCOMBC.rk.2[[l]]))/
      length(unique(c(ANCOMBC.rk.1[[l]], ANCOMBC.rk.2[[l]])))
    ANCOMBC.pr[[l]] <- percent.tx
  }
  for(l in 1:length(replicates)){
    Jaccard <- c(Jaccard, ANCOMBC.pr[[l]])
    methods <- c(methods, "ANCOM-BC")
  }
  study <- c(rep("loso", 5), "all")
  ############################### ALDEx2 ####################################
  Aldex2.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_order/Models_original/Aldex2.model.all.Rdata")
      Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
      Aldex2.rk.1[[k]] <- which(Aldex2.fdr <= fdr)
    }else{
      load(paste0(data.loc.2, "Models_original/Aldex2.model.", as.character(l),".Rdata"))
      Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
      Aldex2.rk.1[[k]] <- which(Aldex2.fdr <= fdr)
    }
    k <- k + 1
  }
  Aldex2.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/Aldex2.model.all.Rdata")
      Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
      Aldex2.rk.2[[k]] <- which((Aldex2.fdr)[filters] <= fdr)
    }else{
      load(paste0(data.loc.1, "Models_original/Aldex2.model.", as.character(l),".Rdata"))
      Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
      Aldex2.rk.2[[k]] <- which(Aldex2.fdr[filters] <= fdr)
    }
    k <- k + 1
  }
  Aldex2.pr <- list()
  for(l in 1:length(replicates)){
    percent.tx <- length(intersect(Aldex2.rk.1[[l]], Aldex2.rk.2[[l]]))/
      length(unique(c(Aldex2.rk.1[[l]], Aldex2.rk.2[[l]])))
    Aldex2.pr[[l]] <- percent.tx
  }
  for(l in 1:length(replicates)){
    Jaccard <- c(Jaccard, Aldex2.pr[[l]])
    methods <- c(methods, "ALDEx2")
  }
  study <- c(rep("loso", 5), "all")
  ############################### BW-prop ###################################
  BW.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_order/Models_original/BW.prop.model.all.Rdata")
    }else{
      load(paste0(data.loc.2, "Models_original/BW.prop.model.", as.character(l),".Rdata"))
    }
    BW.rk.1[[k]] <- which(BW.prop.model$q.val <= fdr)
    k <- k + 1
  }

  BW.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/BW.prop.model.all.Rdata")
    }else{
      load(paste0(data.loc.1, "Models_original/BW.prop.model.", as.character(l),".Rdata"))
    }
    BW.rk.2[[k]] <- which(BW.prop.model$q.val[filters] <= fdr)
    k <- k + 1
  }
  BW.pr <- list()
  for(l in 1:length(replicates)){
    percent.tx <- length(intersect(BW.rk.1[[l]], BW.rk.2[[l]]))/
      length(unique(c(BW.rk.1[[l]], BW.rk.2[[l]])))
    BW.pr[[l]] <- percent.tx
  }
  for(l in 1:length(replicates)){
    Jaccard <- c(Jaccard, BW.pr[[l]])
    methods <- c(methods, "BW")
  }
  study <- c(rep("loso", 5), "all")
  ############################### CLR-LASSO #################################
  clrlasso.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_order/Models_original/CLR.lasso.model.all.Rdata")
    }else{
      load(paste0(data.loc.2, "Models_original/CLR.lasso.model.", as.character(l), ".Rdata"))
    }
    clrlasso.rk.1[[k]] <- names(which(clrlasso$glmnet.fit$beta[,clrlasso$lambda == clrlasso$lambda.min]!=0))
    k <- k + 1
  }
  clrlasso.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/CLR.lasso.model.all.Rdata")
    }else{
      load(paste0(data.loc.1,"Models_original/CLR.lasso.model.", as.character(l), ".Rdata"))
    }
    clrlasso.tmp <-  clrlasso$glmnet.fit$beta[filters,]
    clrlasso.rk.2[[k]] <- names(which(clrlasso.tmp[,clrlasso$lambda == clrlasso$lambda.min] != 0))
    k <- k + 1
  }
  clrlasso.pr <- list()
  for(l in 1:length(replicates)){
    percent.tx <- length(intersect(clrlasso.rk.1[[l]], clrlasso.rk.2[[l]]))/
      length(unique(c(clrlasso.rk.1[[l]], clrlasso.rk.2[[l]])))
    clrlasso.pr[[l]] <- percent.tx
  }
  for(l in 1:length(replicates)){
    Jaccard <- c(Jaccard, clrlasso.pr[[l]])
    methods <- c(methods, "CLR-LASSO")
  }
  study <- c(rep("loso", 5), "all")
  df1 <- data.frame(Jaccard = Jaccard, methods = methods, study = study, type = "Original")

  # Batch-corrected data
  ############################### ANCOM-BC ##################################
  Jaccard <- NULL
  Ham_dist <- NULL
  methods <- NULL
  study <- NULL
  ANCOMBC.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_order/Models_batch_corrected/ANCOMBC.model.all.Rdata")
      ANCOMBC.rk.1[[k]] <- which(ANCOMBC.model$res$q_val$labels1 <= fdr)
    }else{
      load(paste0(data.loc.2, "Models_batch_corrected/ANCOMBC.model.", as.character(l),".Rdata"))
      ANCOMBC.rk.1[[k]] <- which(ANCOMBC.model$res$q_val$labels1 <= fdr)
    }
    k <- k + 1
  }
  ANCOMBC.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_batch_corrected/ANCOMBC.model.all.Rdata")
      ANCOMBC.rk.2[[k]] <- which(ANCOMBC.model$res$q_val$labels1[filters] <= fdr)
    }else{
      load(paste0(data.loc.1, "Models_batch_corrected/ANCOMBC.model.", as.character(l),".Rdata"))
      ANCOMBC.rk.2[[k]] <- which(ANCOMBC.model$res$q_val$labels1[filters] <= fdr)
    }
    k <- k + 1
  }
  ANCOMBC.pr <- list()
  for(l in 1:length(replicates)){
    percent.tx <- length(intersect(ANCOMBC.rk.1[[l]], ANCOMBC.rk.2[[l]]))/
      length(unique(c(ANCOMBC.rk.1[[l]], ANCOMBC.rk.2[[l]])))
    ANCOMBC.pr[[l]] <- percent.tx
  }
  for(l in 1:length(replicates)){
    Jaccard <- c(Jaccard, ANCOMBC.pr[[l]])
    methods <- c(methods, "ANCOM-BC")
  }
  study <- c(rep("loso", 5), "all")
  ############################### ALDEx2 ###################################
  Aldex2.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_order/Models_batch_corrected/Aldex2.model.all.Rdata")
      Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
      Aldex2.rk.1[[k]] <- which(Aldex2.fdr <= fdr)
    }else{
      load(paste0(data.loc.2, "Models_batch_corrected/Aldex2.model.", as.character(l),".Rdata"))
      Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
      Aldex2.rk.1[[k]] <- which(Aldex2.fdr <= fdr)
    }
    k <- k + 1
  }
  Aldex2.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_batch_corrected/Aldex2.model.all.Rdata")
      Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
      Aldex2.rk.2[[k]] <- which((Aldex2.fdr)[filters] <= fdr)
    }else{
      load(paste0(data.loc.1, "Models_batch_corrected/Aldex2.model.", as.character(l),".Rdata"))
      Aldex2.fdr <- Aldex2.model$glmfit$`labels1:pval.BH`
      Aldex2.rk.2[[k]] <- which(Aldex2.fdr[filters] <= fdr)
    }
    k <- k + 1
  }
  Aldex2.pr <- list()
  for(l in 1:length(replicates)){
    percent.tx <- length(intersect(Aldex2.rk.1[[l]], Aldex2.rk.2[[l]]))/
      length(unique(c(Aldex2.rk.1[[l]], Aldex2.rk.2[[l]])))
    Aldex2.pr[[l]] <- percent.tx
  }
  for(l in 1:length(replicates)){
    Jaccard <- c(Jaccard, Aldex2.pr[[l]])
    methods <- c(methods, "ALDEx2")
  }
  study <- c(rep("loso", 5), "all")
  ################################ BW-prop ##################################
  BW.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_order/Models_batch_corrected/BW.prop.model.all.Rdata")
    }else{
      load(paste0(data.loc.2, "Models_batch_corrected/BW.prop.model.", as.character(l),".Rdata"))
    }
    BW.rk.1[[k]] <- which(BW.prop.model$q.val <= fdr)
    k <- k + 1
  }
  BW.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_batch_corrected/BW.prop.model.all.Rdata")
    }else{
      load(paste0(data.loc.1, "Models_batch_corrected/BW.prop.model.", as.character(l),".Rdata"))
    }
    BW.rk.2[[k]] <- which(BW.prop.model$q.val[filters] <= fdr)
    k <- k + 1
  }
  BW.pr <- list()
  for(l in 1:length(replicates)){
    percent.tx <- length(intersect(BW.rk.1[[l]], BW.rk.2[[l]]))/
      length(unique(c(BW.rk.1[[l]], BW.rk.2[[l]])))
    BW.pr[[l]] <- percent.tx
  }
  for(l in 1:length(replicates)){
    Jaccard <- c(Jaccard, BW.pr[[l]])
    methods <- c(methods, "BW")
  }
  study <- c(rep("loso", 5), "all")
  ############################### CLR-LASSO #################################
  clrlasso.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_order/Models_batch_corrected/CLR.lasso.model.all.Rdata")
    }else{
      load(paste0(data.loc.2, "Models_batch_corrected/CLR.lasso.model.", as.character(l), ".Rdata"))
    }
    clrlasso.rk.1[[k]] <- names(which(clrlasso$glmnet.fit$beta[,clrlasso$lambda == clrlasso$lambda.min]!=0))
    k <- k + 1
  }
  clrlasso.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_batch_corrected/CLR.lasso.model.all.Rdata")
    }else{
      load(paste0(data.loc.1,"Models_batch_corrected/CLR.lasso.model.", as.character(l), ".Rdata"))
    }
    clrlasso.tmp <-  clrlasso$glmnet.fit$beta[filters,]
    clrlasso.rk.2[[k]] <- names(which(clrlasso.tmp[,clrlasso$lambda == clrlasso$lambda.min] != 0))
    k <- k + 1
  }
  clrlasso.pr <- list()
  for(l in 1:length(replicates)){
    percent.tx <- length(intersect(clrlasso.rk.1[[l]], clrlasso.rk.2[[l]]))/
      length(unique(c(clrlasso.rk.1[[l]], clrlasso.rk.2[[l]])))
    clrlasso.pr[[l]] <- percent.tx
  }
  for(l in 1:length(replicates)){
    Jaccard <- c(Jaccard, clrlasso.pr[[l]])
    methods <- c(methods, "CLR-LASSO")
  }
  study <- c(rep("loso", 5), "all")
  df2 <- data.frame(Jaccard = Jaccard, methods = methods, study = study, type = "Batch-corrected")

  # Prepare plotting
  df <- rbind(df1, df2)
  df$methods <- factor(df$methods,  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"))
  df$type <- factor(df$type, c("Original", "Batch-corrected"))

  Jac_maxs <- NULL
  Jac_mins <- NULL
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
  pdf("./figures/figure2_e3.pdf",         # File name
      width = 4.72, height = 4.45, # Width and height in inches
      bg = "white")

  df  %>% ggplot(aes(x=methods, y=Jaccard, group = type, color = methods, shape = study)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    scale_shape_manual(breaks = c("loso","all"), values=c(20, 23)) + ylim(0,1) +
    geom_errorbar( aes(x = methods,
                       ymin = Jac_min,
                       ymax = Jac_max,
                       color= methods,
                       linetype = type),
                   width = 0,
                   position = position_dodge(0.5)) +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" )
    ) + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
              panel.background = element_rect(fill = 'white'),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",  text = element_text(size=15))
  ## pdf 5.20 * 5.20

  # Closing the graphical device
  dev.off()
