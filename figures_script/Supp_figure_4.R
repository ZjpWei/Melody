# =============================================== #
#             Supplementary figure 4              #
# =============================================== #

  # Packages ----
  library("tidyverse")
  library("corrplot")
  library("glmnet")
  
  # General---
  load("./CRC_meta_analysis/data/tax.Rdata")
  rownames(tax) <- tax$OTUID
  tax$mOTU <- gsub("[][]","",tax$mOTU)
  Taxa_lst <- tax
  
  # Load Models
  ## Melody list ----
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/Melody.model.all.Rdata")
  taxa.id <- gsub("[][]", "",unlist(regmatches(names(Melody.model$coef),
                                               gregexpr("\\[.*?\\]",names(Melody.model$coef)))))
  Taxa_lst$Melody <- NA
  Taxa_lst[taxa.id,"Melody"] <- rank(abs(Melody.model$coef)) * sign(Melody.model$coef)
  
  ## ALDEx2 list ----
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/Aldex2.model.all.Rdata")
  Taxa_lst$ALDEx2 <- NA
  taxa.id <- gsub("[][]", "",unlist(regmatches(rownames(Aldex2.model$glmfit),
                                               gregexpr("\\[.*?\\]",rownames(Aldex2.model$glmfit)))))
  Taxa_lst[taxa.id,"ALDEx2"] <- rank(-Aldex2.model$glmfit[, "labels1:pval"]) * sign(Aldex2.model$glmfit[, "labels1:Est"])
  
  # Top selected taxa
  ## ANCOMBC list ----
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/ANCOMBC.model.all.Rdata")
  Taxa_lst$"ANCOM-BC" <- NA
  taxa.id <- gsub("[][]", "",unlist(regmatches(ANCOMBC.model$res$q_val$taxon,
                                               gregexpr("\\[.*?\\]",ANCOMBC.model$res$q_val$taxon))))
  Taxa_lst[taxa.id,"ANCOM-BC"] <- rank(-ANCOMBC.model$res$p_val$labels1) * sign(ANCOMBC.model$res$W$labels1)
  
  ## BW list ----
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/BW.prop.model.all.Rdata")
  Taxa_lst$BW <- NA
  taxa.id <- gsub("[][]", "",unlist(regmatches(rownames(BW.prop.model),
                                               gregexpr("\\[.*?\\]",rownames(BW.prop.model)))))
  Taxa_lst[taxa.id,"BW"] <- rank(-BW.prop.model[, "p.val"]) *sign(BW.prop.model[, "statis"])
  
  ## CLR-LASSO list ----
  load("./CRC_meta_analysis/CRC_all_K401/Models_original/CLR.lasso.model.all.Rdata")
  lentmp <- c()
  for(l in 1:ncol(clrlasso$glmnet.fit$beta)){
    lentmp <- c(lentmp, sum(clrlasso$glmnet.fit$beta[,l]!=0))
  }
  tax.id <- gsub("[][]", "",unlist(regmatches(rownames(clrlasso$glmnet.fit$beta),
                                              gregexpr("\\[.*?\\]",rownames(clrlasso$glmnet.fit$beta)))))
  Taxa_lst$`CLR-LASSO` <- NA
  Taxa_lst[taxa.id,"CLR-LASSO"]<- rank(abs(clrlasso$glmnet.fit$beta[, min(which(lentmp >= 51))])) *
    sign(clrlasso$glmnet.fit$beta[, min(which(lentmp >= 51))])
  
  ## Heatmap ----
  species.heatmap.orderd <- names(which(Melody.model$coef!=0))[order(Melody.model$coef[names(which(Melody.model$coef!=0))],decreasing = T)]
  tmp.id <- gsub("[][]", "",unlist(regmatches(species.heatmap.orderd,gregexpr("\\[.*?\\]",species.heatmap.orderd))))
  
  data <- as.matrix(Taxa_lst[tmp.id, c(9:13)])
  data[abs(data) < 350] <- NA 
  rownames(data) <- species.heatmap.orderd
  
  sign.data <- sign(data)
  data <- (abs(data) - 349) * sign.data
  
  # Rename the taxa name
  renames <- sub("^\\S+\\s+",
                 '', tax[gsub("[][]", "",unlist(regmatches(rownames(data) ,
                                                           gregexpr("\\[.*?\\]",rownames(data) )))),"mOTU"])
  # 1. "sp.." to "species" and remove following part
  renames <- gsub(" sp..*"," species",renames)
  renames <- gsub(" subsp..*", " subspecies", renames)
  renames <- gsub(" 3_1_57FAA_CT1", "", renames)
  # 2. remove [C] and remove []
  renames <- gsub("C","",renames)
  renames <- gsub("\\[|\\]","",renames)
  # 3. "gen." to "genus"
  renames <- gsub("gen..*","genus",renames)
  # 4. shorten Fusobacterium
  renames <- gsub("Fusobacterium","F.",renames)
  # Upper unknown
  renames <- str_trim(sub("unknown", "Unknown", renames))
  renames <- paste0(renames, " (",sub("]",")", sub(".*v2_", "", species.heatmap.orderd),")"))
  renames <- factor(renames, levels=renames)
  rownames(data) <- renames

  pdf("./figures/Supp_figure_4.pdf",         # File name
      width = 20, height = 10, # Width and height in inches
      bg = "white")
  
  corrplot((t(data) / apply(abs(data), 2, max, na.rm = TRUE)), 
           type="full" , tl.col= "black", tl.srt=90, tl.cex = 1,
           na.label = "square",na.label.col = "white",
           col = c("#4393C7","#D66051"))
  
  dev.off()
  
