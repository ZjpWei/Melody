# =============================================== #
#       Figure 5 Stability: All v.s. LOSO         #
# =============================================== #

  # Packages ----
  library("ggplot2")
  library("tidyverse")
  library("glmnet")
  library("ggtext")

  # General ----
  rm(list = ls())
  data.all.loc <- "./CRC_meta_analysis/CRC_all_K401/"
  data.loso.loc <- "./CRC_meta_analysis/CRC_loso_K401/"
  load(paste0(data.all.loc, "prepare_data/data.rel.all.Rdata"))
  L <- length(data.rel)
  taxa.names <- colnames(data.rel[[1]]$Y)
  fdr <- 0.01

  # All v.s. LOSO ----
  ## Original data ----
  ### Melody ----
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
  ### ANCOM-BC ----
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
  ### ALDEx2 ----
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
  ### BW-prop ----
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
  ### CLR-LASSO ----
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

  ## Batch-corrected data ----
  ### ANCOM-BC ----
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
  ### ALDEx2 ----
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
  ### BW-prop ----
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
  ### CLR-LASSO ----
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

  ## Prepare plotting ----
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

  p1 <- df  %>% ggplot(aes(x=methods, y=Jaccard, group = type, color = methods)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) + ylim(0,1) +
    geom_errorbar(
      aes(x = methods,ymin = Jac_min, ymax = Jac_max, color= methods,linetype = type),
      width = 0, position = position_dodge(0.5)) +
    ggtitle("All vs LOSO data") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" )
    ) +  theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               plot.title = element_text(hjust = 0.5),
               axis.title = element_blank(),
               axis.text.y = element_text(size = 15),
               axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
               panel.background = element_rect(fill = 'white'),
               axis.line = element_line(colour = "black"),
               panel.border = element_rect(colour = "black", fill=NA),
               legend.position = "bottom",  text = element_text(size=15),
               legend.box="vertical",
               legend.title = element_text(size = 16),
               legend.text = element_text(size = 16),
               legend.background = element_rect(fill = 'NA'),
               legend.key =  element_rect(colour = "transparent", fill = "transparent"),
               legend.key.size = unit(0.5, 'cm'), #change legend key size
               legend.key.height = unit(0.5, 'cm'), #change legend key height
               legend.key.width = unit(1, 'cm'),
               strip.text = element_text(size = 16)) +
    labs(linetype="Pooled data")

# =============================================== #
#        Figure 5 Stability: K401 v.s. K849       #
# =============================================== #
  
  # K401 v.s. K849 ----
  rm(list = setdiff(ls(), "p1") )
  load("./CRC_meta_analysis/data/tax.Rdata")
  load("./CRC_meta_analysis/CRC_all_K849/prepare_data/data.rel.all.Rdata")
  taxa.name1 <- colnames(data.rel[[1]]$Y)
  load("./CRC_meta_analysis/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  taxa.name2 <- colnames(data.rel[[1]]$Y)
  filters <- match(taxa.name2, taxa.name1)
  tax.species <- gsub("[][]", "",unlist(regmatches(colnames(data.rel[[1]]$Y),
                                                   gregexpr("\\[.*?\\]",colnames(data.rel[[1]]$Y)))))
  tax.sub <- tax[match(tax.species, tax$OTUID),]
  fdr <- 0.01
  replicates <- c(1:length(data.rel),"all")
  data.loc.1 <- "./CRC_meta_analysis/CRC_loso_K849/"
  data.loc.2 <- "./CRC_meta_analysis/CRC_loso_K401/"
  
  ## Original data ----
  ### Melody ----
  Jaccard <- NULL
  Ham_dist <- NULL
  methods <- NULL
  study <- NULL
  meta.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/Melody.model.all.Rdata")
      tmpid.1 <- names(Melody.model$coef)
    }else{
      load(paste0(data.loc.2, "Models_original/Melody.model.", as.character(l),".Rdata"))
      tmpid.1 <- names(Melody.model$coef)
    }
    tmp.id <- Melody.model$coef != 0
    meta.rk.1[[k]] <- names(which(tmp.id))
    k <- k + 1
  }
  meta.rk.2 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K849/Models_original/Melody.model.all.Rdata")
      tmpid.2 <- names(Melody.model$coef)
    }else{
      load(paste0(data.loc.1, "Models_original/Melody.model.", as.character(l),".Rdata"))
      tmpid.2 <- names(Melody.model$coef)
    }
    meta.rk.2[[k]] <- names(which(Melody.model$coef[intersect(tmpid.1, tmpid.2)] != 0))
    k <- k + 1
  }
  meta.pr <- list()
  meta.hm <- list()
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
  ### ANCOM-BC ----
  ANCOMBC.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/ANCOMBC.model.all.Rdata")
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
      load("./CRC_meta_analysis/CRC_all_K849/Models_original/ANCOMBC.model.all.Rdata")
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
  ### ALDEx2 ----
  Aldex2.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/Aldex2.model.all.Rdata")
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
      load("./CRC_meta_analysis/CRC_all_K849/Models_original/Aldex2.model.all.Rdata")
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
  ### BW-prop ----
  BW.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/BW.prop.model.all.Rdata")
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
      load("./CRC_meta_analysis/CRC_all_K849/Models_original/BW.prop.model.all.Rdata")
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
  ### CLR-LASSO ----
  clrlasso.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_original/CLR.lasso.model.all.Rdata")
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
      load("./CRC_meta_analysis/CRC_all_K849/Models_original/CLR.lasso.model.all.Rdata")
    }else{
      load(paste0(data.loc.1,"Models_original/CLR.lasso.model.", as.character(l), ".Rdata"))
    }
    clrlasso.tmp <- clrlasso$glmnet.fit$beta[filters,]
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
  
  ## Batch-corrected data ----
  ### ANCOM-BC ----
  Jaccard <- NULL
  Ham_dist <- NULL
  methods <- NULL
  study <- NULL
  ANCOMBC.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_batch_corrected/ANCOMBC.model.all.Rdata")
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
      load("./CRC_meta_analysis/CRC_all_K849/Models_batch_corrected/ANCOMBC.model.all.Rdata")
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
  ### ALDEx2 ----
  Aldex2.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_batch_corrected/Aldex2.model.all.Rdata")
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
      load("./CRC_meta_analysis/CRC_all_K849/Models_batch_corrected/Aldex2.model.all.Rdata")
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
  ### BW-prop ----
  BW.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_batch_corrected/BW.prop.model.all.Rdata")
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
      load("./CRC_meta_analysis/CRC_all_K849/Models_batch_corrected/BW.prop.model.all.Rdata")
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
  ### CLR-LASSO ----
  clrlasso.rk.1 <- list()
  k <- 1
  for(l in replicates){
    if(l == "all"){
      load("./CRC_meta_analysis/CRC_all_K401/Models_batch_corrected/CLR.lasso.model.all.Rdata")
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
      load("./CRC_meta_analysis/CRC_all_K849/Models_batch_corrected/CLR.lasso.model.all.Rdata")
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
  
  ## Prepare plotting ----
  df <- rbind(df1, df2)
  df$methods <- factor(df$methods, c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"))
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
  
  p2 <- df  %>% ggplot(aes(x=methods, y=Jaccard, group = type, color = methods, shape = study)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    scale_shape_manual(breaks = c("loso","all"), values=c(20, 23)) + ylim(0,1) +
    geom_errorbar( aes(x = methods,ymin = Jac_min, ymax = Jac_max, color= methods,
                       linetype = type), width = 0, position = position_dodge(0.5)) +
    ggtitle("<i>K</i> = 401 vs <i>K</i> = 849 species") +
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" )
    ) + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_markdown(hjust = 0.5),
              axis.title = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
              panel.background = element_rect(fill = 'white'),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",  text = element_text(size=15))

# =============================================== #
#     Figure 5 Stability: K401 v.s. K267(order)   #
# =============================================== #
  
  # K401 v.s. K267 ----
  rm(list = setdiff(ls(), c("p1", "p2")) )
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
  
  ## Original data ----
  ### Melody ----
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
  ### ANCOM-BC ----
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
  ### ALDEx2 ----
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
  ### BW-prop ----
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
  ### CLR-LASSO ----
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
  
  ## Batch-corrected data ----
  ### ANCOM-BC ----
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
  ### ALDEx2 ----
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
  ### BW-prop ----
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
  ### CLR-LASSO ----
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
  
  ## Prepare plotting ----
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
  
  p3 <- df  %>% ggplot(aes(x=methods, y=Jaccard, group = type, color = methods, shape = study)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    scale_shape_manual(breaks = c("loso","all"), values=c(20, 23)) + ylim(0,1) +
    
    geom_errorbar( aes(x = methods,
                       ymin = Jac_min,
                       ymax = Jac_max,
                       color= methods,
                       linetype = type),
                   width = 0,
                   position = position_dodge(0.5)) +
    ggtitle("<i>K</i> = 401 vs <i>K</i> = 267 species")+
    scale_color_manual(
      breaks =  c("Melody", "ALDEx2", "ANCOM-BC","BW","CLR-LASSO"),
      values = c("red", "green" ,"blue","grey", "orange" )
    ) + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              plot.title = element_markdown(hjust = 0.5),
              axis.text.y = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 16),
              panel.background = element_rect(fill = 'white'),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.position = "none",  text = element_text(size=15))
  
  # Generate figures ----
  pdf("./figures/figure5.pdf", width = 14, height = 6, bg = "white")
  
  ggpubr::annotate_figure(ggpubr::ggarrange(p1, p2, p3, nrow = 1, ncol = 3, 
                                            common.legend = TRUE, legend = "bottom"), 
                          left = grid::textGrob("Jaccard Index", rot = 90, vjust = 0.5, hjust = 0,
                                                gp = grid::gpar(cex = 1.5)))
  
  dev.off()
  
