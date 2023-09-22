  ##############################################
#                                                #
#       Get top X species for each methods       #
#                 (Random Split)                 #
#                                                #
  ##############################################
  
  rm(list = ls())
  # Packages
  library("glmnet")
  
  data.loc <- "./Prediction/Random_split/"
  if(!dir.exists(paste0(data.loc, "index"))){
    dir.create(paste0(data.loc, "index"))
  }
  if(!dir.exists(paste0(data.loc, "index_batch"))){
    dir.create(paste0(data.loc, "index_batch"))
  }
  load(paste0(data.loc, "prepare_data/data.rel.all.Rdata"))
  L <- length(data.rel)
  taxa.num <- c(1:16) * 5
  taxa.name <- colnames(data.rel[[1]]$Y)
  
  # Original data
  for(s in 1:100){
    load(paste0(data.loc, "Models/Melody.model.",as.character(s),".Rdata"))
    index <- list()
    for(ll in 1:length(taxa.num)){
      index[[ll]] <- names(which(Melody.model$coef[, taxa.num[ll]] != 0))
    }
    names(index) <- as.character(taxa.num)
    save(index, file = paste0(data.loc, "index/index_Melody." , as.character(s),".Rdata"))
    
    load(paste0(data.loc, "Models/Aldex2.model.",as.character(s),".Rdata"))
    Aldex2.names <- rownames(Aldex2.model$glmfit)
    Aldex2.rank <- rank(Aldex2.model$glmfit$`labels1:pval`)
    index <- list()
    for(ll in 1:length(taxa.num)){
      index[[ll]] <- Aldex2.names[Aldex2.rank <= taxa.num[ll]]
    }
    names(index) <- as.character(taxa.num)
    save(index, file = paste0(data.loc, "index/index_Aldex2.", as.character(s) ,".Rdata"))
    
    load(paste0(data.loc, "Models/ANCOMBC.model.",as.character(s),".Rdata"))
    ANCOMBC.rank <- rank(ANCOMBC.model$res$p_val$labels1) 
    index <- list()
    for(ll in 1:length(taxa.num)){
      index[[ll]] <- ANCOMBC.model$res$p_val$taxon[ANCOMBC.rank <= taxa.num[ll]]
    }
    names(index) <- as.character(taxa.num)
    save(index, file = paste0(data.loc, "index/index_ANCOM.",as.character(s) ,".Rdata"))
    
    load(paste0(data.loc, "Models/BW.prop.model.",as.character(s),".Rdata"))
    index <- list()
    BW.rank <- rank(BW.prop.model$p.val)
    for(ll in 1:length(taxa.num)){
      index[[ll]] <- taxa.name[BW.rank <= taxa.num[ll]]
    }
    names(index) <- as.character(taxa.num)
    save(index, file = paste0(data.loc, "index/index_BW.", as.character(s) ,".Rdata"))
    
    load(paste0(data.loc, "Models/CLR.lasso.model.",as.character(s),".Rdata"))
    lentmp <- c()
    for(l in 1:ncol(clrlasso$glmnet.fit$beta)){
      lentmp <- c(lentmp, sum(clrlasso$glmnet.fit$beta[,l]!=0))
    }
    index <- list()
    for(ll in 1:length(taxa.num)){
      tmp.id <- rank(-abs(clrlasso$glmnet.fit$beta[,min(which(lentmp >= taxa.num[ll]))])) <= taxa.num[ll]
      index[[ll]] <- taxa.name[tmp.id]
    }
    names(index) <- as.character(taxa.num)
    save(index, file = paste0(data.loc, "index/index_clrlasso.", as.character(s) ,".Rdata"))
  }
  
  # Batch-corrected data
  for(s in 1:100){
    load(paste0(data.loc, "Models_batch/Aldex2.model.",as.character(s),".Rdata"))
    Aldex2.names <- rownames(Aldex2.model$glmfit)
    Aldex2.rank <- rank(Aldex2.model$glmfit$`labels1:pval`)
    index <- list()
    for(ll in 1:length(taxa.num)){
      index[[ll]] <- Aldex2.names[Aldex2.rank <= taxa.num[ll]]
    }
    save(index, file = paste0(data.loc, "index_batch/index_Aldex2.", as.character(s) ,".Rdata"))
    
    load(paste0(data.loc, "Models_batch/ANCOMBC.model.",as.character(s),".Rdata"))
    ANCOMBC.rank <- rank(ANCOMBC.model$res$p_val$labels1) 
    index <- list()
    for(ll in 1:length(taxa.num)){
      index[[ll]] <- ANCOMBC.model$res$p_val$taxon[ANCOMBC.rank <= taxa.num[ll]]
    }
    names(index) <- as.character(taxa.num)
    save(index, file = paste0(data.loc, "index_batch/index_ANCOM.",as.character(s) ,".Rdata"))
    
    load(paste0(data.loc, "Models_batch/BW.prop.model.",as.character(s),".Rdata"))
    index <- list()
    BW.rank <- rank(BW.prop.model$p.val)
    for(ll in 1:length(taxa.num)){
      index[[ll]] <- taxa.name[BW.rank <= taxa.num[ll]]
    }
    names(index) <- as.character(taxa.num)
    save(index, file = paste0(data.loc, "index_batch/index_BW.", as.character(s) ,".Rdata"))
    
    load(paste0(data.loc, "Models_batch/CLR.lasso.model.",as.character(s),".Rdata"))
    lentmp <- c()
    for(l in 1:ncol(clrlasso$glmnet.fit$beta)){
      lentmp <- c(lentmp, sum(clrlasso$glmnet.fit$beta[,l]!=0))
    }
    index <- list()
    for(ll in 1:length(taxa.num)){
      tmp.id <- rank(-abs(clrlasso$glmnet.fit$beta[,min(which(lentmp >= taxa.num[ll]))])) <= taxa.num[ll]
      index[[ll]] <- taxa.name[tmp.id]
    }
    names(index) <- as.character(taxa.num)
    save(index, file = paste0(data.loc, "index_batch/index_clrlasso.", as.character(s) ,".Rdata"))
  }
  
