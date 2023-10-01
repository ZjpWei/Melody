  ##############################################
#                                                #
#         Figure 2 (c) method overlapping        #
#                                                #
  ##############################################

  # Packages
  library("pheatmap")

  rm(list = ls())
  # Load models
  load("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.ref1.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.ref2.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.ref3.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.ref4.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.ref5.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.mixref.Rdata")

  # Prepare plotting
  index <- list()
  index[[1]] <- names(which(Melody.model.ref1$coef!=0))
  index[[2]] <- names(which(Melody.model.ref2$coef!=0))
  index[[3]] <- names(which(Melody.model.ref3$coef!=0))
  index[[4]] <- names(which(Melody.model.ref4$coef!=0))
  index[[5]] <- names(which(Melody.model.ref5$coef!=0))
  index[[6]] <- names(which(Melody.model.mixref$coef!=0))
  Jaccard <- matrix(0,6,6)
  for(l in 1:6){
    for(ll in 1:6)
      Jaccard[l,ll] <- length(intersect(index[[l]], index[[ll]]))/length(unique(c(index[[l]], index[[ll]])))
  }
  rownames(Jaccard) <- paste0(c("ref1", "ref2", "ref3","ref4", "ref5", "ref.mix"),
                              " (", unlist(lapply(index, length)),")")

  pdf("./figures/figure2_c_sensitivity.pdf",         # File name
      width = 8, height = 5.68, # Width and height in inches
      bg = "white")

  # Creating a plot
  pheatmap::pheatmap(round(Jaccard,2), cluster_rows = FALSE, cluster_cols = FALSE,
                     breaks = c(1:100)/100, display_numbers = TRUE,
                     fontsize_number = 22, fontsize_row=20)


  # Closing the graphical device
  dev.off()
