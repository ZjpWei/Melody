# =============================================== #
#             Supplementary figure 10             #
# =============================================== #

  # Packages ----
  library("RColorBrewer")
  library("stringr")
  library("tidyverse")
  library("ggtext")
  library("cowplot")
  library("grid")
  library("gridExtra")
  library("ggplot2")
  library("miMeta")

  ## Read data
  rm(list = ls())
  load("./Metabolites/Processed_data/processed_data.Rdata")

  ## process data
  cluster <- lapply(metadata, function(d){
    d <- d %>% dplyr::transmute(Sample, Subject) %>% data.frame()
    rownames(d) <- NULL
    d <- d %>% tibble::column_to_rownames("Sample")
    return(d)
  })

  Metabolite.ids <- unique(common.pairs$Compound)
  Genera.id <- unique(common.pairs$Taxon)

  covariates.disease <- lapply(metadata[datasets], function(d){
    case.names <- c("Healthy", "Control", "Normal", "H", "0", "nonIBD", "NA")
    return(d %>% tibble::rownames_to_column(var = "IDS") %>%
             dplyr::transmute(Sample, Study.Group =  case_when(Study.Group %in% case.names ~ "0",
                                                               .default = "1"))
    )
  })

  #################################### Choose original compounds ####################################
  datasets.cmpd <- c()
  covariates_adjust_lst <- list()
  for(d in datasets){
    if(length(data.for.lm[[d]]$cmpd.name) != length(unique(data.for.lm[[d]]$cmpd.name))){
      otu_data <- matrix(0, nrow = nrow(data.for.lm[[d]]$rel.abd), ncol = length(Genera.id),
                         dimnames = list(rownames(data.for.lm[[d]]$rel.abd), Genera.id))
      datasets.cmpd <- c(datasets.cmpd, d)
      cmpd.name <- names(which(table(data.for.lm[[d]]$cmpd.name) > 1))
      cmpd.dat <- data.for.lm[[d]]$rel.cmpd[,data.for.lm[[d]]$cmpd.name %in% cmpd.name]
      otu_data[,colnames(data.for.lm[[d]]$rel.abd)] <- data.for.lm[[d]]$rel.abd %>% as.matrix()
      sample.kep <- intersect(rownames(otu_data), rownames(cmpd.dat))

      cmpd.lst <- c()
      for(l in cmpd.name){
        cmpd.lst <- rbind(cmpd.lst,
                          cbind(colnames(data.for.lm[[d]]$rel.cmpd)[data.for.lm[[d]]$cmpd.name == l],
                                l))
      }
      colnames(cmpd.lst) <- c("Orig.Compound", "Compound")
      cmpd.lst <- data.frame(cmpd.lst)
      covariates_adjust_lst[[d]] <- covariates.disease[[d]][match(sample.kep,covariates.disease[[d]]$Sample),] %>%
        data.frame(row.names = NULL) %>% tibble::column_to_rownames(var = "Sample")
    }
  }

  for(d in datasets){
    if(d %in% datasets.cmpd){
      load(paste0("./Metabolites/Metabolite_select/Compound_", d, ".Rdata"))
      remove.cmpd <- setdiff(cmpd_sty_scan$Orig.Compound, cmpd_sty_scan %>%
                               dplyr::filter(selected_num == max(selected_num), .by = Compound) %>%
                               dplyr::filter(row_number() == 1, .by = Compound) %>%
                               pull(Orig.Compound))
      keep.cmpd <- !(colnames(data.for.lm[[d]]$rel.cmpd) %in% remove.cmpd)
      cmpd.dat <- data.for.lm[[d]]$rel.cmpd[,keep.cmpd]
      data.for.lm[[d]]$cmpd.name <- data.for.lm[[d]]$cmpd.name[keep.cmpd]
      colnames(cmpd.dat) <- data.for.lm[[d]]$cmpd.name
      data.for.lm[[d]]$rel.cmpd <- cmpd.dat
    }else{
      colnames(data.for.lm[[d]]$rel.cmpd) <- data.for.lm[[d]]$cmpd.name
    }
  }

  ## Reform data
  otu_data_lst <- list()
  cmpd_data_lst <- list()
  cluster_data_lst <- list()
  covariates_adjust_lst <- list()
  for(d in datasets){
    if("HMDB0001161" %in% colnames(data.for.lm[[d]]$rel.cmpd)){
      sample.kep <- intersect(rownames(data.for.lm[[d]]$rel.abd), rownames(data.for.lm[[d]]$rel.cmpd))
      covariates_adjust_lst[[d]] <- covariates.disease[[d]][match(sample.kep, covariates.disease[[d]]$Sample),] %>%
        data.frame(row.names = NULL) %>%
        tibble::column_to_rownames(var = "Sample")
      otu_data_lst[[d]] <- data.for.lm[[d]]$rel.abd[sample.kep,]
      cmpd_data_lst[[d]] <- data.for.lm[[d]]$rel.cmpd[sample.kep, ] %>% dplyr::transmute(HMDB0001161 = HMDB0001161)
      cluster_data_lst[[d]] <- (cluster[[d]])[sample.kep,]
      names(cluster_data_lst[[d]]) <- sample.kep
    }
  }
  covariates_adjust_lst$POYET_BIO_ML_2019 <- NULL

  # Melody heatmap
  p_melody <- list()
  rm.gn <- c()
  for(k in 1:2){
    otu_data_lst_sub <- list()
    cmpd_data_lst_sub <- list()
    for(d in names(otu_data_lst)){
      otu_data_lst_sub[[d]] <- (otu_data_lst[[d]])[,!(colnames(otu_data_lst[[d]]) %in% rm.gn)]
      cmpd_data_lst_sub[[d]] <- cmpd_data_lst[[d]] %>% dplyr::transmute(HMDB0001161 = HMDB0001161)
    }

    null.obj <- melody.null.model(rel.abd = otu_data_lst_sub, covariate.adjust = covariates_adjust_lst)

    summary.stat.study <- melody.get.summary(null.obj = null.obj,
                                             covariate.interest = cmpd_data_lst_sub,
                                             cluster = cluster_data_lst)

    Melody.model <- melody.meta.summary(summary.stats = summary.stat.study,
                                        tune.path = "sequence",
                                        output.best.one = TRUE)

    len.melody <- sum(Melody.model$HMDB0001161$coef!=0)

    rm.gn <- names(which(rank(-Melody.model$HMDB0001161$coef) <= 5))

    # Main ----
    L <- length(summary.stat.study)
    species.name <- names(Melody.model$HMDB0001161$coef)
    tab.AA <- NULL
    for(d in 1:L){
      tmp.AA <- rep(NA, length(species.name))
      names(tmp.AA) <- species.name
      tmp.AA[names(summary.stat.study[[d]]$est[,1])] <- summary.stat.study[[d]]$est[,1] + Melody.model$HMDB0001161$delta[d]
      tmp.AA[names(summary.stat.study[[d]]$est[,1])] <- tmp.AA[names(summary.stat.study[[d]]$est[,1])]
      tab.AA <- rbind(tab.AA, tmp.AA)
    }
    est.vals <- data.frame(all = Melody.model$HMDB0001161$coef)
    rownames(tab.AA) <- paste0("MTBL", as.character(match(names(otu_data_lst), datasets)))
    colnames(tab.AA) <- species.name
    tab.AA <- as.data.frame(t(tab.AA))
    ref.studies <- paste0("MTBL", as.character(match(names(otu_data_lst), datasets)))
    species.heatmap <- rownames(est.vals)[which(est.vals$all != 0)]
    est.vals.signed <- est.vals[species.heatmap,"all", drop=FALSE]
    species.heatmap.orderd <- c(rownames(est.vals.signed[order(est.vals.signed[,"all"],decreasing = T),,drop=FALSE]))
    renames <- species.heatmap.orderd

    # Prepare plotting
    est.vals.plot <- est.vals[species.heatmap.orderd,"all"]
    tab.AA.plot <- tab.AA[species.heatmap.orderd,]
    mx <- 1
    brs <- seq(-mx, mx, by=1e-3)
    num.col.steps <- length(brs)-1
    n <- floor(0.5*num.col.steps)
    col.hm <- c(rev(colorRampPalette(brewer.pal(9, 'Blues'))(n)),
                rep('#FFFFFF', num.col.steps-2*n),
                colorRampPalette(brewer.pal(9, 'Reds'))(n))
    tab.AA.plot[tab.AA.plot >= mx] <- mx
    tab.AA.plot[tab.AA.plot <= -mx] <- -mx

    # Heatmap plot ----
    rownames(tab.AA.plot) <- gsub(".*;g__","", rownames(tab.AA.plot))
    renames <- gsub(".*;g__", "", renames)
    plot.single.study.heatmap.ref <- function(x){
      df.plot <- NULL
      for(l in x){
        df.plot <- rbind(df.plot,
                         tibble(species=factor(renames, levels=renames),
                                study = factor(rep(l, nrow(tab.AA.plot)), levels = l),
                                AA=tab.AA.plot[,l]))
      }
      g1 <- df.plot %>%
        ggplot(aes(x=species, y=study, fill=AA)) +
        geom_tile() + theme_minimal() + ylab("Study") +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(size = 15),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect(fill=NULL, colour='black'),
              plot.margin = unit(c(0, 0, 0, 0), 'cm'),
              axis.text.y = element_text(size = 13),
              axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1, size = 12, face = "italic"),
              legend.title = element_text(hjust = 0.5, size = 16),
              legend.text = element_text(size = 15),
              legend.position = "bottom") +
        scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx),
                             breaks = c(-mx,-0.5,0,0.5,mx),
                             labels = as.character(c(-mx,-0.5,0,0.5,mx))) +
        guides(fill = guide_colourbar(barwidth = 12, barheight = 0.6,
                                      title = "Study-specific AA coefficient",
                                      title.position = "top"))
    }
    g1 <- plot.single.study.heatmap.ref(rev(ref.studies))


    # Scatter plot ----
    g2 <- tibble(species=factor(renames, levels=renames),
                 est.vals= est.vals.plot,
                 variable = "AA beta_hat est",
                 colour=est.vals != 0) %>%
      ggplot() + geom_point(aes(x= species, y= est.vals)) +
      theme_classic() + ylab('Meta <br> coefficient') +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 15),
            axis.title.y = element_markdown(size = 15),
            axis.title.x = element_blank(),
            legend.position = "right") +
      scale_y_continuous(limits=c(min(est.vals.plot) - 0.05, max(est.vals.plot) + 0.05), expand = c(0, 0)) +
      scale_x_discrete(position='bottom') +
      scale_fill_manual(values=c('lightgrey', 'darkgrey'), guide="none") +
      geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed")

    # Generate figures ----
    p_melody[[k]] <- plot_grid(g2, g1, nrow=2, align = 'v', rel_heights = c(.25,.75))
  }

  # MMUPHin heatmap
  source("./utility/mmuphin.R")
  load("./Metabolites/Results/Melody.null.model.Rdata")
  otu_data_lst <- NULL
  cmpd_data_lst <- NULL
  cluster_data_lst <- NULL
  covariates_adjust_lst <- NULL
  ref.studies <- c()
  for(d in datasets){
    if("HMDB0001161" %in% colnames(data.for.lm[[d]]$rel.cmpd)){
      sample.kep <- intersect(rownames(data.for.lm[[d]]$rel.abd),
                              rownames(data.for.lm[[d]]$rel.cmpd)[!is.na(data.for.lm[[d]]$rel.cmpd[,"HMDB0001161"])])
      covariates_adjust_lst <- rbind(covariates_adjust_lst, covariates.disease[[d]][match(sample.kep, covariates.disease[[d]]$Sample),] %>%
                                       data.frame(row.names = NULL) %>%
                                       tibble::column_to_rownames(var = "Sample"))

      tmp.otu <- matrix(0, nrow = length(sample.kep), ncol = length(Genera.id),
                        dimnames = list(sample.kep, Genera.id))
      inter.taxa <- colnames(null.obj[[d]]$res)
      tmp.otu[sample.kep, inter.taxa] <- as.matrix(data.for.lm[[d]]$rel.abd[sample.kep,inter.taxa])
      otu_data_lst <- rbind(otu_data_lst, tmp.otu)
      tmp.cmpd <- matrix(0, nrow = length(sample.kep), ncol = 1,
                         dimnames = list(sample.kep, "HMDB0001161"))
      tmp.cmpd[sample.kep, "HMDB0001161"] <- as.matrix(data.for.lm[[d]]$rel.cmpd[sample.kep,"HMDB0001161"])
      cmpd_data_lst <- rbind(cmpd_data_lst, tmp.cmpd)
      cluster_data_lst <- rbind(cluster_data_lst, cluster[[d]])
      ref.studies <- c(ref.studies, rep(d, nrow(tmp.otu)))
    }
  }
  cluster_data_lst <- cluster_data_lst[rownames(otu_data_lst),]
  meta.data <- cbind(cluster_data_lst, covariates_adjust_lst, ref.studies, cmpd_data_lst)
  colnames(meta.data) <- c("Subject", "Study.Group", "Study" ,"HMDB0001161")
  meta.data$Study <- factor(meta.data$Study)

  p_mmuphin <- list()
  rm.gm <- c()
  for(k in 1:2){
    if(length(unique(meta.data$Subject)) == length(meta.data$Subject)){
      MMUPHin.model <- fit_metaAnalysis(feature.abd = t(round(otu_data_lst[,!(colnames(otu_data_lst) %in% rm.gm)])),
                                        data = meta.data,
                                        test_variable = "HMDB0001161",
                                        batch_variable = "Study",
                                        covariates = "Study.Group",
                                        covariates.random = NULL,
                                        moderator_variables = NULL)
    }else{
      MMUPHin.model <- fit_metaAnalysis(feature.abd = t(round(otu_data_lst[,!(colnames(otu_data_lst) %in% rm.gm)])),
                                        data = meta.data,
                                        test_variable = "HMDB0001161",
                                        batch_variable = "Study",
                                        covariates = "Study.Group",
                                        covariates.random = "Subject",
                                        moderator_variables = NULL)
    }
    MMUPHin.model$maaslin_fits <- MMUPHin.model$maaslin_fits[order(match(names(MMUPHin.model$maaslin_fits), datasets))]

    summary.stat.study <- list()
    for(l in 1:length(MMUPHin.model$maaslin_fits)){
      summary.stat.study[[l]] <- list(est = MMUPHin.model$maaslin_fits[[l]]$coef)
      names(summary.stat.study[[l]]$est) <- MMUPHin.model$maaslin_fits[[l]]$feature
    }
    MMUPHin.coef <- MMUPHin.model$meta_fits$coef
    names(MMUPHin.coef) <- MMUPHin.model$meta_fits$feature
    MMUPHin.coef[p.adjust(MMUPHin.model$meta_fits$pval, method = "fdr") > 0.05] <- 0
    len.melody <- sum(MMUPHin.coef!=0)
    rm.gm <- names(which(rank(-MMUPHin.coef) <= 5))

    # Main ----
    L <- length(summary.stat.study)
    species.name <- gsub(".*;g__","g__", names(MMUPHin.coef))
    tab.AA <- NULL
    for(l in 1:L){
      tmp.AA <- rep(NA, length(species.name))
      names(tmp.AA) <- species.name
      tmp.AA[names(summary.stat.study[[l]]$est)] <- summary.stat.study[[l]]$est
      tab.AA <- rbind(tab.AA, tmp.AA)
    }
    est.vals <- data.frame(all = MMUPHin.coef)

    rownames(tab.AA) <- paste0("MTBL", as.character(sort(match(levels(meta.data$Study), datasets))))
    ref.studies <- paste0("MTBL", as.character(sort(match(levels(meta.data$Study), datasets))))
    tab.AA <- as.data.frame(t(tab.AA))
    species.heatmap <- rownames(est.vals)[which(est.vals$all != 0)]
    est.vals.signed <- est.vals[species.heatmap,"all", drop=FALSE]
    species.heatmap.orderd <- c(rownames(est.vals.signed[order(est.vals.signed[,"all"],decreasing = T),,drop=FALSE]))
    renames <- gsub(".*;g__","g__",species.heatmap.orderd)

    # Prepare plotting
    est.vals.plot <- est.vals[species.heatmap.orderd,"all"]
    tab.AA.plot <- tab.AA[species.heatmap.orderd,]
    mx <- 0.02
    brs <- seq(-mx, mx, by=1e-4)
    num.col.steps <- length(brs)-1
    n <- floor(0.5*num.col.steps)
    col.hm <-   col.hm <- c(rev(colorRampPalette(brewer.pal(9, 'Purples'))(n)),
                            rep('#FFFFFF', num.col.steps-2*n),
                            colorRampPalette(c("#F9F9FC", "#FFE4E4", "#FFC5C5",
                                               "#FFA6A6", "#FF8787", "#FF6868",
                                               "#FF4848", "#FF2929", "#FF0000"))(n))

    tab.AA.plot[tab.AA.plot >= mx] <- mx
    tab.AA.plot[tab.AA.plot <= -mx] <- -mx

    # Heatmap plot ----
    rownames(tab.AA.plot) <- gsub(".*;g__","", rownames(tab.AA.plot))
    renames <- gsub("g__", "", renames)
    plot.single.study.heatmap.ref <- function(x){
      df.plot <- NULL
      for(l in x){
        df.plot <- rbind(df.plot,
                         tibble(species=factor(renames, levels=renames),
                                study = factor(rep(l, nrow(tab.AA.plot)), levels = l),
                                AA=tab.AA.plot[,l]))
      }
      g1 <- df.plot %>%
        ggplot(aes(x=species, y=study, fill=AA)) +
        geom_tile() + theme_minimal() + ylab("Study") +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(size = 15),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect(fill=NULL, colour='black'),
              plot.margin = unit(c(0, 0, 0, 0), 'cm'),
              axis.text.y = element_text(size = 13),
              axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1, size = 12, face = "italic"),
              legend.title = element_text(hjust = 0.5, size = 16),
              legend.text = element_text(size = 15),
              legend.position = "bottom") +
        scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx),
                             breaks = c(-mx,-mx/2,0,mx/2,mx),
                             labels = as.character(c(-mx,-mx/2,0,mx/2,mx))) +
        guides(fill = guide_colourbar(barwidth = 12, barheight = 0.6,
                                      title = "Study-specific coefficient",
                                      title.position = "top"))
    }
    g1 <- plot.single.study.heatmap.ref(rev(ref.studies))

    # Scatter plot ----
    g2 <- tibble(species=factor(renames, levels=renames),
                 est.vals= est.vals.plot,
                 variable = "AA beta_hat est",
                 colour=est.vals != 0) %>%
      ggplot() + geom_point(aes(x= species, y= est.vals)) +
      theme_classic() + ylab('Meta <br> coefficient') +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 15),
            axis.title.y = element_markdown(size = 15),
            axis.title.x = element_blank(),
            legend.position = "right") +
      scale_y_continuous(limits=c(min(est.vals.plot) - 0.005, max(est.vals.plot) + 0.005), expand = c(0, 0)) +
      scale_x_discrete(position='bottom') +
      scale_fill_manual(values=c('lightgrey', 'darkgrey'), guide="none") +
      geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed")

    # Generate figures ----
    p_mmuphin[[k]] <- plot_grid(g2, g1, nrow=2, align = 'v', rel_heights = c(.25,.75))
  }


  pdf("./figures/FigS10_melody_bf.pdf", width = 6.75, height = 5.10, bg = "white")

  p_melody[[1]]

  dev.off()

  pdf("./figures/FigS10_melody_af.pdf", width = 6.75, height = 5.10, bg = "white")

  p_melody[[2]]

  dev.off()


  pdf("./figures/FigS10_mmuphin_bf.pdf", width = 6.75, height = 5.10, bg = "white")

  p_mmuphin[[1]]

  dev.off()

  pdf("./figures/FigS10_mmphin_af.pdf", width = 6.75, height = 5.10, bg = "white")

  p_mmuphin[[2]]

  dev.off()


