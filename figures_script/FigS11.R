# =============================================== #
#             Supplementary figure 11             #
# =============================================== #

  # Packages ----
  library("tidyverse")
  library("ggplot2")

  # Main ----
  rm(list = ls())
  load("./Metabolites/Results/Melody.Jaccard.Rdata")
  load("./Metabolites/Results/MMUPHin.Jaccard.Rdata")
  rm.cp.lst <- colnames(MMUPHin.af.coef)[-1]

  ## Jaccard plot
  len.compare <- c()
  Jaccard_m1 <- c()
  for(l in rm.cp.lst){
    names1 <- setdiff(Melody.bf$feature[Melody.bf[[l]]!=0], Melody.bf$feature[rank(-Melody.bf[[l]]) <= 5])
    names2 <- Melody.af.best$feature[Melody.af.best[[l]]!=0 & !is.na(Melody.af.best[[l]])]
    len.compare <- rbind(len.compare, c(length(names1), length(names2)))
    Jaccard_m1 <- c(Jaccard_m1, length(intersect(names1,names2))/length(unique(c(names1,names2))))
  }
  rownames(len.compare) <- rm.cp.lst
  names(Jaccard_m1) <- rm.cp.lst


  len.compare <- c()
  Jaccard_v1 <- c()
  for(l in rm.cp.lst){
    bf.tmp.coef <- MMUPHin.bf.coef[[l]]
    bf.tmp.coef[p.adjust(MMUPHin.bf.pval[[l]], method = "BH") > 0.05] <- 0
    af.tmp.coef <- MMUPHin.af.coef[[l]]
    af.tmp.coef[p.adjust(MMUPHin.af.pval[[l]], method = "BH") > 0.05] <- 0

    names1 <- setdiff(MMUPHin.bf.coef$feature[bf.tmp.coef!=0], MMUPHin.bf.coef$feature[rank(-bf.tmp.coef) <= 5])
    names2 <- MMUPHin.af.coef$feature[af.tmp.coef!=0 & !is.na(af.tmp.coef)]
    len.compare <- rbind(len.compare, c(length(names1), length(names2)))
    Jaccard_v1 <- c(Jaccard_v1, length(intersect(names1,names2))/length(unique(c(names1,names2))))
  }
  rownames(len.compare) <- rm.cp.lst
  names(Jaccard_v1) <- rm.cp.lst


  df_scatter <- data.frame(MMUPHin = Jaccard_v1) %>% tibble::rownames_to_column("Compounds") %>%
    dplyr::left_join( data.frame(Melody = Jaccard_m1) %>% tibble::rownames_to_column("Compounds"), by = "Compounds")

  pdf("./figures/FigS11.pdf", width = 6.27, height = 5.94, bg = "white")

  df_scatter %>% ggplot(aes(x = MMUPHin, y = Melody)) +

    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("MMUPHin Jaccard index") + ylab("Melody Jaccard index") +
    geom_point() + theme_minimal() + ylim(0,1) + xlim(0,1) +
    theme(axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17),
          title = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 0.1, size = 18),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          legend.title = element_text(hjust = 0.5, size = 15),
          legend.text = element_text(size = 15),
          legend.position = "none")

  dev.off()


