# =============================================== #
#         Figure 3 (b): Method overlapping        #
# =============================================== #

  # Packages ----
  library("pheatmap")
  library("glmnet")

  # General ----
  rm(list = ls())
  load("./CRC/Processed_data/data.org.K401.Rdata")
  L <- length(data.rel)

  # Load models and data
  load("./CRC/Models_org/Melody.model.Rdata")
  load("./CRC/Models_org/BW.prop.model.Rdata")
  load("./CRC/Models_org/ANCOMBC2.model.Rdata")
  load("./CRC/Models_org/Aldex2.model.Rdata")
  load("./CRC/Models_org/MMUPHin.model.Rdata")
  load("./CRC/Models_org/CLR.lasso.model.Rdata")

  # Get selected taxa for each method (other methods match the number selected by Melody)
  index <- list()
  taxa.id <- gsub("[][]", "",unlist(regmatches(names(which(Melody.model$disease$coef!=0)),
                             gregexpr("\\[.*?\\]",names(which(Melody.model$disease$coef!=0))))))
  len.tax <- length(taxa.id)
  index[[1]] <- taxa.id

  aldex2.top <- rownames(Aldex2.model$glmfit)[rank(Aldex2.model$glmfit[,c("labels1:pval")]) <= len.tax]
  taxa.id <- gsub("[][]", "",unlist(regmatches(aldex2.top, gregexpr("\\[.*?\\]",aldex2.top))))
  index[[2]] <- taxa.id

  ancom.top <- ANCOMBC2.model$res$taxon[rank(ANCOMBC2.model$res$p_labels1) <= len.tax]
  taxa.id <- gsub("[][]", "",unlist(regmatches(ancom.top, gregexpr("\\[.*?\\]",ancom.top))))
  index[[3]] <- taxa.id

  mmuphin.top <- MMUPHin.model$meta_fits$feature[rank(MMUPHin.model$meta_fits$pval) <= len.tax]
  taxa.id <- gsub("[][]", "",unlist(regmatches(mmuphin.top, gregexpr("\\[.*?\\]",mmuphin.top))))
  index[[4]] <- taxa.id

  BW.top <- rownames(BW.prop.model)[rank(BW.prop.model$p.val) <= len.tax]
  taxa.id <- gsub("[][]", "",unlist(regmatches(BW.top, gregexpr("\\[.*?\\]",BW.top))))
  index[[5]] <- taxa.id

  len.tmp <- c()
  for(ll in 1:ncol(clrlasso$glmnet.fit$beta)){
    len.tmp <- c(len.tmp, sum(clrlasso$glmnet.fit$beta[,ll]!=0))
  }
  taxa.id <- names(which(rank(-abs(clrlasso$glmnet.fit$beta[,which(len.tmp >= len.tax)[1]])) <= len.tax))
  taxa.id <- gsub("[][]", "",unlist(regmatches(taxa.id, gregexpr("\\[.*?\\]",taxa.id))))
  index[[6]] <- taxa.id

  # prepare plotting ----
  Heatmap <- matrix(0,6,6)
  for(l in 1:length(index)){
    for(ll in 1:length(index)) {
      Heatmap[l, ll] <- length(intersect(index[[l]], index[[ll]]))/
        length(unique(c(index[[l]], index[[ll]])))
    }
  }
  rownames(Heatmap) <- c("Melody", "ALDEx2", "ANCOM-BC2", "MMUPHin", "BW", "CLR-LASSO")

  # Generate figures ----
  pdf("./figures/Fig3_b.pdf", width = 8, height = 5.68, bg = "white")

  pheatmap::pheatmap(round(Heatmap,2), cluster_rows = TRUE, cluster_cols = TRUE,
                     color = (colorRampPalette(c("#2A788EFF", "#7AD151FF", "#FDE725FF"))(100))[20:100],
                     breaks = c(20:100)/100, display_numbers = TRUE,fontsize_number = 22, fontsize = 18)

  dev.off()
