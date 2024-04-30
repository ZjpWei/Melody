# =============================================== #
#             Supplementary figure 9              #
# =============================================== #

  # Packages ----
  library(plotly)
  library(tidyr)
  library(ggplot2)
  library(cluster)
  library(ggdendro)
  library(grid)
  library(randomcoloR)
  library(ggtext)
  
  rm(list = ls())
  load("./Metabolites/Processed_data/processed_data.Rdata")
  load("./Metabolites/Results/Melody.meta.all.Rdata")
  load("./Metabolites/Results/MMUPHin.meta.all.Rdata")
  load("./Metabolites/Results/summary.stat.study.all.Rdata")
  feature.ids <- rownames(Melody.meta.model$HMDB0000008$coef)
  Melody_coef_mat <- matrix(0, nrow = length(feature.ids), ncol = length(Melody.meta.model),
                            dimnames = list(feature.ids, names(Melody.meta.model)))
  for(ids in names(Melody.meta.model)){
    best.coef <- Melody.meta.model[[ids]]$coef[,which.min(Melody.meta.model[[ids]]$dev + Melody.meta.model[[ids]]$ic)]
    Melody_coef_mat[rownames(Melody.meta.model[[ids]]$coef), ids] <- best.coef
  }
  Melody_coef_mat[is.na(Melody_coef_mat)] <- 0
  MMUPHin_coef_mat <- MMUPHin.meta.coef %>% tibble::column_to_rownames(var = "feature")
  MMUPHin_pval_mat <- MMUPHin.meta.pval %>% tibble::column_to_rownames(var = "feature")
  MMUPHin_qval_mat <- apply(MMUPHin_pval_mat, 2, p.adjust, method = "fdr")
  MMUPHin_coef_mat[MMUPHin_qval_mat > 0.05] <- 0
  MMUPHin_coef_mat[is.na(MMUPHin_coef_mat)] <- 0
  
  sample.id <- as.vector(unlist(lapply(data.for.lm, function(d){rownames(d$rel.abd)})))
  taxa <- NULL
  for(d in names(data.for.lm)){
    taxa <- c(taxa, colnames(data.for.lm[[d]]$rel.abd))
  }
  taxa <- unique(taxa)
  Otu_big_mat <- matrix(0, nrow = length(sample.id), ncol = length(taxa),
                        dimnames = list(sample.id, taxa))
  for(d in names(data.for.lm)){
    Otu_big_mat[rownames(data.for.lm[[d]]$rel.abd),colnames(data.for.lm[[d]]$rel.abd)] <- as.matrix(data.for.lm[[d]]$rel.abd)
  }
  tax.prevalence <- colMeans(Otu_big_mat / rowSums(Otu_big_mat))
  names(tax.prevalence) <- gsub(".*;g__","g__", names(tax.prevalence))
  
  df_genera <- data.frame(Melody = rowSums(Melody_coef_mat!=0)) %>% 
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(MMUPHin = rowSums(MMUPHin_coef_mat!=0)) %>% 
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>% 
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>% 
    dplyr::arrange(desc(Melody)) %>% 
    dplyr::mutate(group = dplyr::case_when(Genera %in%  c("g__Ruminococcus_B",
                                                          "g__Faecalimonas",
                                                          "g__Sutterella",         
                                                          "g__Clostridium_Q",
                                                          "g__Anaerostipes",
                                                          "g__Phocaeicola",
                                                          "g__Roseburia",
                                                          "g__Erysipelatoclostridium",
                                                          "g__Enterocloster",
                                                          "g__Blautia",
                                                          "g__Clostridium",
                                                          "g__Veillonella",
                                                          "g__Haemophilus_D") ~ "unique genera",
    TRUE ~ "common genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "unique genera" ~ gsub("g__","", Genera),
                                           TRUE ~ "")) %>% 
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) 
  
  df_cmpd <- data.frame(Melody = colSums(Melody_coef_mat != 0)) %>% 
    tibble::rownames_to_column(var = "Compounds") %>% dplyr::left_join(
      data.frame(MMUPHin = colSums(MMUPHin_coef_mat != 0)) %>% 
        tibble::rownames_to_column(var = "Compounds"), by = "Compounds"
    ) %>% dplyr::arrange(desc(Melody)) %>% 
    dplyr::mutate(HMDB.Class = "TRUE")
  
  df_genera$hjust[df_genera$label == "Blautia"] <- 0.8
  df_genera$hjust[df_genera$label == "Phocaeicola"] <- 0.2
  df_genera$vjust[df_genera$label == "Faecalimonas"] <- -3
  df_genera$hjust[df_genera$label == "Faecalimonas"] <- 0.5
  df_genera$vjust[df_genera$label == "Veillonella"] <- -1.5
  df_genera$hjust[df_genera$label == "Veillonella"] <- 0
  df_genera$vjust[df_genera$label == "Sutterella"] <- -1
  df_genera$vjust[df_genera$label == "Haemophilus_D"] <- -2
  df_genera$vjust[df_genera$label == "Clostridium"] <- -1
  df_genera$vjust[df_genera$label == "Erysipelatoclostridium"] <- -3.5
  df_genera$vjust[df_genera$label == "Clostridium_Q"] <- 1.4
  df_genera$hjust[df_genera$label == "Clostridium_Q"] <- 0
  
  p_genus <- df_genera %>% ggplot(aes(x = MMUPHin, y = Melody, color = group, label = label)) + 
    geom_point(aes(size = prevalence)) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in MMUPHin results")+ ylab("Count in Melody results") +
    theme_minimal() + ylim(0, 200) + xlim(0, 200) + ggtitle("") +
    scale_color_manual(values=c("gray", "darkolivegreen3")) + 
    scale_size_manual(values = c(4,3,2,1), 
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_genera$hjust, vjust=df_genera$vjust, size = 3, color = "black") + 
    theme(axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17), 
          title = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 0.1, size = 18),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = c(0.85,0.2)) + labs(size="Trimmed average <br> proportion") +
    guides(color = "none")
  
 
  p_cmpd <- df_cmpd %>% ggplot(aes(x = MMUPHin, y = Melody, color = HMDB.Class)) + 
    geom_point() +  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in MMUPHin results")+ 
    ylab("Count in Melody results") +
    theme_minimal() + ylim(0, 70) + xlim(0, 70) + ggtitle("") +
    scale_color_manual(values="orange") + 
    theme(axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17), 
          title = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 0.1, size = 18),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_text(hjust = 0.5, size = 15),
          legend.text = element_text(size = 15),
          legend.position = "none")

  pdf("./figures/FigS9_genus.pdf", width = 7.04, height = 6.64, bg = "white")
  
  p_genus
  
  dev.off()
  
  pdf("./figures/FigS9_cmpd.pdf", width = 7.04, height = 6.04, bg = "white")
  
  p_cmpd
  
  dev.off()
  
  