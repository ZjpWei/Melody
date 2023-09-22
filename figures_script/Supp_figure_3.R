  ##############################################
#                                                #
#             Supplementary figure 3             #
#                                                #
  ##############################################
  
  library("ggplot2")
  
  rm(list = ls())
  # General
  load("./CRC_Real/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  L <- length(data.rel)
  # Train Melody model with sequence:
  # Run models
  L <- length(data.rel)
  Study <- NULL
  Group <- NULL
  rel.abd <- NULL
  for(l in 1:length(data.rel)){
    Study <- c(Study, rep(as.character(l), nrow(data.rel[[l]]$Y)))
    Group <- c(Group, as.character(data.rel[[l]]$X))
    rel.abd <- rbind(rel.abd, data.rel[[l]]$Y)
  }
  rel.abd <- t(rel.abd)
  meta <- data.frame(Sample_ID = colnames(rel.abd), Study = Study, Group = Group)
  rownames(meta) <- meta$Sample_ID
  
  # generate summary statistics
  summary.stat.study <- melody.get.summary(rel.abd = rel.abd,
                                           sample.data = meta,
                                           sample.id = "Sample_ID",
                                           study = "Study",
                                           disease = "Group",
                                           ref = "Coprococcus catus [ref_mOTU_v2_4874]")
  
  # meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study, tune.type = "HBIC",
                                      tune.path = "sequence", tune.size.sequence = 1:50,
                                      ouput.best.one = FALSE)
  
  load("./data/tax.Rdata")
  load("./CRC_Real/CRC_all_K401/Models/BW.prop.model.all.Rdata")
  load("./CRC_Real/CRC_all_K401/Models/ANCOMBC.model.all.Rdata")
  load("./CRC_Real/CRC_all_K401/Models/Aldex2.model.all.Rdata")
  load("./CRC_Real/CRC_all_K401/Models/CLR.lasso.model.all.Rdata")
  
  # Prepare plotting
  for(tops in c(20, 30, 40, 50)){
    index <- list()
    taxa.id <- gsub("[][]", "",unlist(regmatches(names(which(Melody.model$coef[,tops]!=0)),
                   gregexpr("\\[.*?\\]",names(which(Melody.model$coef[,tops]!=0))))))
    len.tax <- length(taxa.id)
    index[[1]] <- taxa.id
    ###
    aldex2.top <- rownames(Aldex2.model$glmfit)[rank(Aldex2.model$glmfit[,c("labels1:pval")]) <= len.tax]
    taxa.id <- gsub("[][]", "",unlist(regmatches(aldex2.top, gregexpr("\\[.*?\\]",aldex2.top))))
    index[[2]] <- taxa.id
    ###
    ancom.top <- ANCOMBC.model$res$p_val$taxon[rank(ANCOMBC.model$res$p_val$labels1) <= len.tax]
    taxa.id <- gsub("[][]", "",unlist(regmatches(ancom.top, gregexpr("\\[.*?\\]",ancom.top))))
    index[[3]] <- taxa.id
    ###
    BW.top <- rownames(BW.prop.model)[rank(BW.prop.model$p.val) <= len.tax]
    taxa.id <- gsub("[][]", "",unlist(regmatches(BW.top, gregexpr("\\[.*?\\]",BW.top))))
    index[[4]] <- taxa.id
    ###
    len.tmp <- c()
    for(ll in 1:ncol(clrlasso$glmnet.fit$beta)){
      len.tmp <- c(len.tmp, sum(clrlasso$glmnet.fit$beta[,ll]!=0))
    }
    taxa.id <- names(which(rank(-abs(clrlasso$glmnet.fit$beta[,which(len.tmp >= len.tax)[1]])) <= len.tax))
    taxa.id <- gsub("[][]", "",unlist(regmatches(taxa.id, gregexpr("\\[.*?\\]",taxa.id))))
    index[[5]] <- taxa.id
    ###############################################
    index <- index[c(5,1,4,2,3)]
    Heatmap <- matrix(0,5,5)
    for(l in 1:length(index)){
      for(ll in 1:length(index)) {
        Heatmap[l, ll] <- length(intersect(index[[l]], index[[ll]]))/length(unique(c(index[[l]], index[[ll]])))
      }
    }
    rownames(Heatmap) <- c("CLR-LASSO", "Melody", "BW", "ALDEx2", "ANCOM-BC")
    pdf(paste0("./figures/Supp_3_top", tops,".pdf"),         # File name
        width = 6, height = 5, # Width and height in inches
        bg = "white")    
    
    # Creating a plot
    pheatmap::pheatmap(Heatmap, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(1:100)/100, 
                       display_numbers = TRUE,fontsize_number = 20, fontsize = 15, legend = FALSE)
    
    
    # Closing the graphical device
    dev.off() 
  }
  