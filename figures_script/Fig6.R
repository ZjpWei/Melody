# =============================================== #
#     Figure 6: Analysis on Metabolites data      #
# =============================================== #

  # Packages ----
  library(plotly)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(cluster)
  library(ggdendro)
  library(grid)
  library(randomcoloR)
  library(ggtext)
  library(RColorBrewer)
  library(cowplot)
  
  # General ----
  rm(list = ls())
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
  
  df1 <- data.frame(Melody = rowSums(Melody_coef_mat!=0)) %>% 
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(MMUPHin = rowSums(MMUPHin_coef_mat!=0)) %>% 
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>% 
    tibble::column_to_rownames("Genera") %>% dplyr::arrange(desc(Melody)) %>% 
    dplyr::slice(which(Melody >= 70))
  
  df2 <- data.frame(Melody = colSums(Melody_coef_mat!=0)) %>% 
    tibble::rownames_to_column(var = "Compounds") %>% dplyr::left_join(
      data.frame(MMUPHin = colSums(MMUPHin_coef_mat!=0)) %>% 
        tibble::rownames_to_column(var = "Compounds"), by = "Compounds"
    ) %>% tibble::column_to_rownames("Compounds") %>% dplyr::arrange(desc(Melody)) %>%
    dplyr::slice(which(Melody >= 30))
  
  df3 <- data.frame(Melody = rowSums(Melody_coef_mat!=0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(MMUPHin = rowSums(MMUPHin_coef_mat!=0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera"
    )  %>% tibble::column_to_rownames("Genera") %>% dplyr::arrange(desc(MMUPHin)) %>% 
    dplyr::slice(which(MMUPHin >= 70))
  
  selected.tax.melody <- rownames(df1)
  selected.cmpd <- rownames(df2)
  selected.tax.mmuphin <- rownames(df3)

  source("./utility/hmdb_utils.R")
  hmdb.ids <- selected.cmpd
  hmdb.data <- get.hmdb.data.by.ids(hmdb.ids)
  hmdb.data <- hmdb.data %>% rename(Compound = HMDB)
  
  ############################# MMUPHin ############################# 
  # CMPD
  cmpd.level <- data.frame(Compound = selected.cmpd) %>% left_join(hmdb.data, by = "Compound") %>%
    mutate(HMDB.Class2 = HMDB.Super.Class) %>%
    mutate(HMDB.Class2 = factor(HMDB.Class2, 
                                levels = c("Benzenoids",            
                                           "Lipids and lipid-like molecules",         
                                           "Organic acids and derivatives",           
                                           "Organic nitrogen compounds",             
                                           "Organic oxygen compounds",                
                                           "Organoheterocyclic compounds",   
                                           "Phenylpropanoids and polyketides"),
                                labels = c("Others", 
                                           "Lipids and lipid-like molecules",  
                                           "Organic acids and derivatives",
                                           "Organic nitrogen compounds",  
                                           "Others",
                                           "Others",
                                           "Others")))
  ## Compute FDR
  Tax.level <- data.frame(Taxon = selected.tax.mmuphin) %>% 
    tidyr::separate(Taxon, into = c("Taxon.Domain","Taxon.Phylum", "Taxon.Class",
                                    "Taxon.Order","Taxon.Family","Taxon.Genus"), 
                    sep = ";", fill = "right", remove = FALSE) %>%
    dplyr::mutate(Taxon.Order.Full = Taxon.Order)
  
  Tax.level$Taxon.Order.Full <- factor(Tax.level$Taxon.Order.Full, 
                                       levels = c("o__Oscillospirales",
                                                  "o__Lachnospirales",
                                                  "o__Bacteroidales",
                                                  "o__Enterobacterales",
                                                  "o__Christensenellales",
                                                  "o__Veillonellales",
                                                  "o__Verrucomicrobiales",
                                                  "o__TANB77",
                                                  "o__Lactobacillales",
                                                  "o__Erysipelotrichales",
                                                  "o__Clostridiales",
                                                  "o__Burkholderiales",
                                                  "o__Actinomycetales"))
  
  selected.tax.mmuphin <-  gsub(".*;g__","g__", selected.tax.mmuphin)
  rownames(MMUPHin_coef_mat) <-  gsub(".*;g__","g__",  rownames(MMUPHin_coef_mat))
  
  ## Hierarchical cluster
  selected.cmpd.order <- c()
  for(sp.cmpd in sort(unique(as.character(cmpd.level$HMDB.Class2)))){
    Compound.id <- cmpd.level$Compound[sp.cmpd == cmpd.level$HMDB.Class2]
    if(length(Compound.id) < 2){
      selected.cmpd.order <- c(selected.cmpd.order, Compound.id)
    }else{
      sub_coef_mat <- Melody_coef_mat[,Compound.id]
      clusters <- hclust(dist(t(sub_coef_mat)))
      selected.cmpd.order <- c(selected.cmpd.order, Compound.id[clusters$order])
    }
  }
  clusters <- hclust(dist(MMUPHin_coef_mat[selected.tax.mmuphin,]))
  selected.tax.order <- selected.tax.mmuphin[clusters$order]
  rownames(Tax.level) <- Tax.level$Taxon.Genus
  Tax.level <- Tax.level[selected.tax.order,]
  rownames(cmpd.level) <- cmpd.level$Compound
  cmpd.level <- cmpd.level[selected.cmpd.order,]
  
  mx <- 0.02
  AA.map <- t(MMUPHin_coef_mat[selected.tax.order, selected.cmpd.order])
  colnames(AA.map) <- gsub("g__", "", colnames(AA.map))
  selected.tax.order <- gsub("g__", "", selected.tax.order)
  Genera_MMUPHin <- selected.tax.order
  
  AA.map[AA.map >= mx] <- mx
  AA.map[AA.map <= -mx] <- -mx
  brs <- seq(-mx, mx, by=1e-6)
  num.col.steps <- length(brs)-1
  n <- floor(0.49*num.col.steps)
  col.hm <- c(rev(colorRampPalette(brewer.pal(9, 'Purples'))(n)),
              rep('#FFFFFF', num.col.steps-2*n),
              colorRampPalette(c("#F9F9FC", "#FFE4E4", "#FFC5C5", 
                                 "#FFA6A6", "#FF8787", "#FF6868", 
                                 "#FF4848", "#FF2929", "#FF0000"))(n))

  plot.single.study.heatmap.mmuphin <- function(x){
    df.plot <- NULL
    for(l in x){
      df.plot <- rbind(df.plot,
                       tibble(species=factor(selected.tax.order, levels=selected.tax.order),
                              cmpd = factor(rep(cmpd.level$HMDB.Name[cmpd.level$Compound == l],
                                                length(selected.tax.order)),
                                            levels = cmpd.level$HMDB.Name[cmpd.level$Compound == l]),
                              AA=AA.map[l,],
                              Taxon.order = Tax.level$Taxon.Order.Full,
                              HMDB.Class2 = rep(as.character(cmpd.level$HMDB.Class2)[l == cmpd.level$Compound],
                                                length(selected.tax.order))))
    }
    
    text.type = rep("italic", length(selected.tax.order))
    text.type[selected.tax.order %in% c("Ruminococcus_B",
                                        "Faecalimonas",
                                        "Sutterella",         
                                        "Clostridium_Q",
                                        "Anaerostipes",
                                        "Phocaeicola",
                                        "Roseburia",
                                        "Erysipelatoclostridium",
                                        "Enterocloster",
                                        "Blautia",
                                        "Clostridium",
                                        "Veillonella",
                                        "Haemophilus_D")] <-"bold.italic"
    
    g1 <- df.plot %>%
      ggplot(aes(x=species, y=cmpd, fill = AA)) +
      geom_tile() + theme_minimal() + ylab("cmpd") +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_rect(fill=NULL, colour='black'),
            plot.margin = unit(c(0, 0, 0, 0), 'cm'),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1, size = 12, face = text.type),
            legend.title = element_text(hjust = 0.5, size = 16),
            legend.text = element_text(size = 15),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vertical",
            strip.text = element_blank()) +
      facet_grid(rows = vars(HMDB.Class2),
                 # col = vars(Taxon.order),
                 scales = "free", space = "free", switch = "y") +
      scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx),
                           breaks = c(-mx,-mx/2,0,mx/2,mx),
                           labels = as.character(c(-mx,-mx/2,0,mx/2,mx))) +
      guides(fill = guide_colourbar(barwidth = 15, barheight = 0.6, title = "Meta coefficient",
                                    title.position = "top"))
  }
  g2 <- plot.single.study.heatmap.mmuphin(rev(selected.cmpd.order))
  
  Tax.level$Taxon.Genus <- gsub("g__", "",   Tax.level$Taxon.Genus)
  Tax.level$Taxon.Genus <- factor(Tax.level$Taxon.Genus, levels =  Tax.level$Taxon.Genus)
  Tax.level$Order <-  Tax.level$Taxon.Order.Full
  Tax.level$Order <- gsub("o__", "",  Tax.level$Order)
  Tax.level$Order[!(Tax.level$Order %in% c("Oscillospirales", "Lachnospirales","Bacteroidales"))] <- "Others"
  
  genera.bar <- ggplot(Tax.level, aes(x = Taxon.Genus, y = 0, fill = Order)) +
    geom_tile()  + theme_minimal() +
    scale_fill_manual(values = c("darkgreen", "orange1","cyan","gray"),
                      breaks = c("Oscillospirales", "Lachnospirales",
                                 "Bacteroidales", "Others")) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none",
          legend.key.size = unit(0.3, 'cm'))  + 
    guides(fill = guide_legend(nrow = 1))
  
  null.bar <- ggplot(Tax.level[1,], aes(x = 0, y = 0, fill = 0)) +
    geom_tile()  + theme_minimal() +
    scale_fill_gradient2(low = "white",mid = "white", high = "white") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='white'),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none")
  
  pdf("./figures/Fig6_MMUPHin.pdf", width = 7.98, height = 13, bg = "white")
  
  plot_grid(plot_grid(null.bar, genera.bar, ncol = 2, align = 'h', rel_widths = c(.3,.7)), g2,
            nrow = 2, align = 'v', rel_heights = c(.015,1))
  
  dev.off()
  
  ############################# Melody ############################# 
  ## Compute FDR
  Tax.level <- data.frame(Taxon = selected.tax.melody) %>% 
    tidyr::separate(Taxon, into = c("Taxon.Domain","Taxon.Phylum", "Taxon.Class",
                                    "Taxon.Order","Taxon.Family","Taxon.Genus"), 
                    sep = ";", fill = "right", remove = FALSE) %>%
    dplyr::mutate(Taxon.Order.Full = Taxon.Order)
  
  Tax.level$Taxon.Order.Full <- factor(Tax.level$Taxon.Order.Full, 
                                       levels = c("o__Oscillospirales",
                                                  "o__Lachnospirales",
                                                  "o__Bacteroidales",
                                                  "o__Enterobacterales",
                                                  "o__Christensenellales",
                                                  "o__Veillonellales",
                                                  "o__Verrucomicrobiales",
                                                  "o__TANB77",
                                                  "o__Lactobacillales",
                                                  "o__Erysipelotrichales",
                                                  "o__Clostridiales",
                                                  "o__Burkholderiales",
                                                  "o__Actinomycetales"))
  
  selected.tax.melody <-  gsub(".*;g__","g__",selected.tax.melody)
  rownames(Melody_coef_mat) <- gsub(".*;g__","g__",rownames(Melody_coef_mat))
  
  ## Hierarchical cluster
  clusters <- hclust(dist(Melody_coef_mat[selected.tax.melody,])) 
  selected.tax.order <- selected.tax.melody[clusters$order]
  rownames(Tax.level) <- Tax.level$Taxon.Genus
  Tax.level <- Tax.level[selected.tax.order,]
  
  mx <- 0.8
  AA.map <- Melody_coef_mat[selected.tax.order, selected.cmpd.order]
  rownames(AA.map) <- gsub("g__", "", rownames(AA.map))
  selected.tax.order <- gsub("g__", "", selected.tax.order)
  Genera_Melody <- selected.tax.order[1:17]
  text.type <- rep("italic", length(selected.tax.order))
  text.type[selected.tax.order %in% setdiff(Genera_Melody, Genera_MMUPHin)] <- "bold.ltalic"
  
  AA.map[AA.map >= mx] <- mx
  AA.map[AA.map <= -mx] <- -mx
  brs <- seq(-mx, mx, by=1e-3)
  num.col.steps <- length(brs)-1
  n <- floor(0.49*num.col.steps)
  col.hm <- c(rev(colorRampPalette(brewer.pal(9, 'Blues'))(n)),
              rep('#FFFFFF', num.col.steps-2*n),
              colorRampPalette(brewer.pal(9, 'Reds'))(n))
  
  
  
  plot.single.study.heatmap.ref <- function(x){
    df.plot <- NULL
    for(l in x){
      df.plot <- rbind(df.plot,
                       tibble(species = factor(selected.tax.order, levels=selected.tax.order),
                              cmpd = factor(rep(cmpd.level$HMDB.Name[cmpd.level$Compound == l],
                                                length(selected.tax.order)),
                                            levels = cmpd.level$HMDB.Name[cmpd.level$Compound == l]),
                              AA=AA.map[,l],
                              # text.type = text.type,
                              Taxon.order = Tax.level$Taxon.Order.Full,
                              HMDB.Class2 = rep(as.character(cmpd.level$HMDB.Class2)[l == cmpd.level$Compound],
                                                length(selected.tax.order))))
    }

    text.type = rep("italic", length(selected.tax.order))
    text.type[selected.tax.order %in% c("Ruminococcus_B",
                                        "Faecalimonas",
                                        "Sutterella",         
                                        "Clostridium_Q",
                                        "Anaerostipes",
                                        "Phocaeicola",
                                        "Roseburia",
                                        "Erysipelatoclostridium",
                                        "Enterocloster",
                                        "Blautia",
                                        "Clostridium",
                                        "Veillonella",
                                        "Haemophilus_D")] <-"bold.italic"
    
    g1 <- df.plot %>%
      ggplot(aes(x=species, y=cmpd, fill = AA)) +
      geom_tile() + theme_minimal() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_rect(fill=NULL, colour='black'),
            plot.margin = unit(c(0, 0, 0, 0), 'cm'),
            axis.text.y = element_text(size = 12),
            axis.text.x =  element_markdown(angle = 90, vjust = 0.5, hjust=1, size = 12, face = text.type),
            legend.title = element_text(hjust = 0.5, size = 16),
            legend.text = element_text(size = 15),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vertical",
            strip.text = element_blank()) +
      facet_grid(rows = vars(HMDB.Class2),
                 scales = "free", space = "free", switch = "y") +
      scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx),
                           breaks = c(-0.8,-0.4,0,0.4,0.8),
                           labels = c("-0.8", "-0.4", "0", "0.4", "0.8")) +
      guides(fill = guide_colourbar(barwidth = 15, barheight = 0.6, title = "Meta coefficient",
                                    title.position = "top"))
  }
  g1 <- plot.single.study.heatmap.ref(rev(selected.cmpd.order))
  
  Tax.level$Taxon.Genus <- gsub("g__", "", Tax.level$Taxon.Genus)
  Tax.level$Taxon.Genus <- factor(Tax.level$Taxon.Genus, levels =  Tax.level$Taxon.Genus)
  Tax.level$Order <- gsub("o__", "", Tax.level$Taxon.Order.Full)
  Tax.level$Order[!(Tax.level$Order %in% c("Oscillospirales", "Lachnospirales","Bacteroidales"))] <- "Others"
  
  genera.bar <- ggplot(Tax.level, aes(x = Taxon.Genus, y = 0, fill = Order)) +
    geom_tile()  + theme_minimal() +
    scale_fill_manual(values = c("darkgreen", "orange1","cyan","gray"),
                      breaks = c("Oscillospirales", "Lachnospirales",
                                 "Bacteroidales", "Others")) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'),
          axis.text.y = element_blank(), 
          axis.text.x = element_blank(),
          legend.position = "none",
          legend.text = element_text(face = "italic"),
          legend.key.size = unit(0.3, 'cm')) + 
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  null.bar <- ggplot(Tax.level[1,], aes(x = 0, y = 0, fill = 0)) +
    geom_tile()  + theme_minimal() +
    scale_fill_gradient2(low = "white",mid = "white", high = "white") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='white'),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'),
          axis.text.y = element_blank(), 
          axis.text.x = element_blank(),
          legend.position = "none")
  
  pdf("./figures/Fig6_Melody.pdf", width = 7.98, height = 13, bg = "white")
  
  plot_grid(plot_grid(null.bar, genera.bar, ncol = 2, align = 'h', rel_widths = c(.3,.7)), g1, 
            nrow = 2, align = 'v', rel_heights = c(.015,1))
  
  dev.off()
  
  plot(clusters, labels = FALSE, hang = -1)
  
