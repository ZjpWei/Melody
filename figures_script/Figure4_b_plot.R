# =============================================== #
#          Figure 4 (b) method overlappinp        #
# =============================================== #

  # Packages ----
  library("pheatmap")
  library("glmnet")

  # General ----
  rm(list = ls())
  load("./CRC_meta_analysis/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  L <- length(data.rel)

  # Load models and data
  load("./CRC_meta_analysis/data/tax.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/Melody.model.all.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/BW.prop.model.all.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/ANCOMBC.model.all.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/Aldex2.model.all.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/CLR.lasso.model.all.Rdata")

  # Get selected taxa for each method (other methods match the number selcted by Melody)
  index <- list()
  taxa.id <- gsub("[][]", "",unlist(regmatches(names(which(Melody.model$coef!=0)),
                             gregexpr("\\[.*?\\]",names(which(Melody.model$coef!=0))))))
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

  # prepare plotting ----
  Heatmap <- matrix(0,5,5)
  for(l in 1:length(index)){
    for(ll in 1:length(index)) {
      Heatmap[l, ll] <- length(intersect(index[[l]], index[[ll]]))/
        length(unique(c(index[[l]], index[[ll]])))
    }
  }
  rownames(Heatmap) <- c("Melody", "ALDEx2", "ANCOM-BC", "BW", "CLR-LASSO")

  # Generate figures ----
  pdf("./figures/figure4_b_comparison.pdf", width = 8, height = 5.68, bg = "white")

  pheatmap::pheatmap(round(Heatmap,2), cluster_rows = TRUE, cluster_cols = TRUE, breaks = c(1:100)/100,
                     display_numbers = TRUE,fontsize_number = 22, fontsize = 18)

  dev.off()
