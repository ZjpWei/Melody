# =============================================== #
#       Figure 3 (c): Reference sensitivity       #
# =============================================== #

  # Packages ----
  library("pheatmap")

  # General ----
  rm(list = ls())
  load("./CRC/sensitivity/Melody.model.ref1.Rdata")
  load("./CRC/sensitivity/Melody.model.ref2.Rdata")
  load("./CRC/sensitivity/Melody.model.ref3.Rdata")
  load("./CRC/sensitivity/Melody.model.ref4.Rdata")
  load("./CRC/sensitivity/Melody.model.ref5.Rdata")
  load("./CRC/sensitivity/Melody.model.mixref.Rdata")

  # Prepare plotting ----
  index <- list()
  index[[1]] <- names(which(Melody.model.ref1$disease$coef!=0))
  index[[2]] <- names(which(Melody.model.ref2$disease$coef!=0))
  index[[3]] <- names(which(Melody.model.ref3$disease$coef!=0))
  index[[4]] <- names(which(Melody.model.ref4$disease$coef!=0))
  index[[5]] <- names(which(Melody.model.ref5$disease$coef!=0))
  index[[6]] <- names(which(Melody.model.mixref$disease$coef!=0))
  Jaccard <- matrix(0,6,6)
  for(l in 1:6){
    for(ll in 1:6)
      Jaccard[l,ll] <- length(intersect(index[[l]], index[[ll]]))/length(unique(c(index[[l]], index[[ll]])))
  }
  rownames(Jaccard) <- paste0(c("ref1", "ref2", "ref3","ref4", "ref5", "ref.mix"),
                              " (", unlist(lapply(index, length)),")")

  # Generate figures ----
  pdf("./figures/Fig3_c.pdf", width = 8, height = 5.68, bg = "white")

  pheatmap::pheatmap(round(Jaccard,2), cluster_rows = FALSE, cluster_cols = FALSE,
                     color = (colorRampPalette(c("#2A788EFF", "#7AD151FF", "#FDE725FF"))(100))[20:100],
                     breaks = c(20:100)/100, display_numbers = TRUE,
                     fontsize_number = 24, fontsize_row=22)

  dev.off()
