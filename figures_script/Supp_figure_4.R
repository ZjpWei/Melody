  ##############################################
#                                                #
#             Supplementary figure 4             #
#                                                #
  ##############################################

  library("RColorBrewer")
  library("stringr")
  library("tidyverse")
  library("ggtext")
  library("cowplot")
  library("grid")
  library("gridExtra")

  # taxa intersection
  rm(list = ls())
  tax.inters <- c()
  for(ss in 1:5){
    load(paste0("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.ref",
                as.character(ss),".Rdata"))
    if(ss == 1){
      Melody.model <- Melody.model.ref1
    }else if(ss == 2){
      Melody.model <- Melody.model.ref2
    }else if(ss == 3){
      Melody.model <- Melody.model.ref3
    }else if(ss ==4){
      Melody.model <- Melody.model.ref4
    }else if(ss == 5){
      Melody.model <- Melody.model.ref5
    }
    tax.inters[[ss]] <- names(which(Melody.model$coef!=0))
  }
  tax.inter <- intersect(tax.inters[[1]], tax.inters[[2]])
  tax.inter <- intersect(tax.inter, tax.inters[[3]])
  tax.inter <- intersect(tax.inter, tax.inters[[4]])
  tax.inter <- intersect(tax.inter, tax.inters[[5]])

  RA.lst <- list()
  AA.lst <- list()
  for(ss in 1:5){
    load(paste0("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.ref",
                as.character(ss),".Rdata"))
    load(paste0("./CRC_meta_analysis/CRC_all_K401/summary.stats/summary.stats.ref",
                as.character(ss),".Rdata"))
    if(ss == 1){
      summary.stat.study <- summary.stat.study.ref1
      Melody.model <- Melody.model.ref1
    }else if(ss == 2){
      summary.stat.study <- summary.stat.study.ref2
      Melody.model <- Melody.model.ref2
    }else if(ss == 3){
      summary.stat.study <- summary.stat.study.ref3
      Melody.model <- Melody.model.ref3
    }else if(ss ==4){
      summary.stat.study <- summary.stat.study.ref4
      Melody.model <- Melody.model.ref4
    }else if(ss == 5){
      summary.stat.study <- summary.stat.study.ref5
      Melody.model <- Melody.model.ref5
    }
    load("./CRC_meta_analysis/data/tax.Rdata")
    len.melody <- length(tax.inter)
    ###############################################
    load("./CRC_meta_analysis/CRC_all_K401/prepare_data/data.rel.all.Rdata")
    L <- length(data.rel)
    species.name <- colnames(data.rel[[1]]$Y)
    ####################
    tab.RA <- NULL
    tab.AA <- NULL
    for(l in 1:L){
      tmp.RA <- rep(NA, ncol(data.rel[[1]]$Y))
      names(tmp.RA) <- colnames(data.rel[[1]]$Y)
      tmp.RA[names(summary.stat.study$summary.stat.study[[l]]$est)] <- summary.stat.study$summary.stat.study[[l]]$est
      tmp.RA[names(Melody.model$delta)[l]] <- 0
      tmp.AA <- tmp.RA + Melody.model$delta[l]
      tab.RA <- rbind(tab.RA, tmp.RA)
      tab.AA <- rbind(tab.AA, tmp.AA)
    }
    est.vals <- data.frame(all = Melody.model$coef)
    rownames(tab.RA) <- paste0("CRC", as.character(1:L))
    colnames(tab.RA) <- species.name
    rownames(tab.AA) <- paste0("CRC", as.character(1:L))
    colnames(tab.AA) <- species.name
    tab.AA <- as.data.frame(t(tab.AA))
    tab.RA <- as.data.frame(t(tab.RA))
    ref.studies <- paste0("CRC", as.character(1:L))
    est.vals.signed <- est.vals[tax.inter,"all", drop=FALSE]
    if(ss == 1){
      tax.inter.orderd <- c(rownames(est.vals.signed[order(est.vals.signed[,"all"],decreasing = T),,drop=FALSE]))
    }

    # prepare plotting
    est.vals.plot <- est.vals[tax.inter.orderd, "all"]
    tab.AA.plot <- tab.AA[tax.inter.orderd,]
    tab.RA.plot <- tab.RA[tax.inter.orderd,]
    mx <- 2.5
    brs = seq(-mx, mx, by=1e-3)
    num.col.steps = length(brs)-1
    n = floor(0.5*num.col.steps)
    col.hm = c(rev(colorRampPalette(brewer.pal(9, 'Blues'))(n)),
               rep('#FFFFFF', num.col.steps-2*n),
               colorRampPalette(brewer.pal(9, 'Reds'))(n))

    tab.AA.plot[tab.AA.plot >= mx] <- mx
    tab.AA.plot[tab.AA.plot <= -mx] <- -mx
    tab.RA.plot[tab.RA.plot >= mx] <- mx
    tab.RA.plot[tab.RA.plot <= -mx] <- -mx
    # Rename the taxa name in plot
    rownames(tax) <- tax$OTUID
    renames <- sub("^\\S+\\s+",
                   '', tax[gsub("[][]", "",unlist(regmatches(tax.inter.orderd,
                                                             gregexpr("\\[.*?\\]",tax.inter.orderd)))),"mOTU"])
    # 1. "sp.." to "sp." and remove following part
    renames <- gsub("sp..*","sp.",renames)
    renames <- sub("subSpecies", "subspecies", renames)
    # 2. remove [C] and remove []
    renames <- gsub("\\[C\\]","",renames)
    renames <- gsub("\\[|\\]","",renames)
    # 3. "gen." to "genus"
    renames <- gsub("gen..*","genus",renames)
    # 4. shorten Fusobacterium
    renames <- gsub("Fusobacterium","F.",renames)
    # Upper unkonwn
    renames <- str_trim(sub("unknown", "Unknown", renames))
    sp_strs <- str_split(renames, " ")
    for(l in 1:length(sp_strs)){
      strs <- sp_strs[[l]]
      for(ll in 1:length(strs)){
        if(strs[[ll]] %in% c("Unknown", "sp.", "genus", "subsp.")){
          strs[[ll]] <- paste0("</i>",strs[[ll]]," <i>")
        }
      }
      strs[[1]] <- paste0(" <i>",strs[[1]])
      strs[[length(strs)]] <- paste0(strs[[length(strs)]],"</i>")
      sp_strs[[l]] <- strs
    }
    renames <- sub("<i></i>", "", unlist(lapply(sp_strs, paste0, collapse=" ")))
    renames <- paste0(renames, " (",sub("]",")", sub(".*v2_", "", tax.inter.orderd),")"))
    renames <- factor(renames, levels=renames)
    ########################################################
    plot.single.study.heatmap <- function(x){
      df.plot <- NULL
      for(l in x){
        df.plot <- rbind(df.plot,
                         tibble(species=factor(renames, levels=renames),
                                study = factor(rep(l, nrow(tab.RA.plot)), levels = l),
                                fc=tab.RA.plot[,l]))
      }
      df.plot$variable <- "RA summary est"
      g1 <- df.plot %>%
        ggplot(aes(x=species, y=study, fill=fc)) +
        geom_tile() + theme_minimal() +
        xlab('Gut microbial species') +
        ylab(paste0("ref", ss)) +
        theme(axis.text.y = element_text(),
              axis.title.y = element_text(size = 15),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect(fill=NULL, colour='black'),
              plot.margin = unit(c(0, 0, 0, 0), 'cm'),
              axis.text = element_blank()) +
        scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx), guide="none")
    }
    ####################################################

    plot.single.study.heatmap.ref <- function(x){
      df.plot <- NULL
      for(l in x){
        df.plot <- rbind(df.plot,
                         tibble(species=factor(renames, levels=renames),
                                study = factor(rep(l, nrow(tab.AA.plot)), levels = l),
                                fc=tab.AA.plot[,l]))
      }
      df.plot$variable <- "AA summary est"
      g1 <- df.plot %>%
        ggplot(aes(x=species, y=study, fill=fc)) +
        geom_tile() + theme_minimal() +
        xlab('Gut microbial species') +
        ylab(paste0("ref",ss)) +
        theme(axis.text.y = element_text(),
              axis.title.y = element_text(size = 15),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect(fill=NULL, colour='black'),
              plot.margin = unit(c(0, 0, 0, 0), 'cm'),
              axis.text = element_blank()) +
        scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx), guide="none")
    }
    # ##############################################################################

    g.lst <- plot.single.study.heatmap(rev(ref.studies))

    g.lst.2 <- plot.single.study.heatmap.ref(rev(ref.studies))

    RA.lst[[ss]] <- g.lst

    AA.lst[[ss]] <- g.lst.2
  }

  # Customizing the output
  pdf("./figures/Supp_4_RA_heatmap.pdf",         # File name
      width = 7.67, height = 3.78, # Width and height in inches
      bg = "white")

  # Creating a plot
  plot_grid(RA.lst[[1]],
            RA.lst[[2]],
            RA.lst[[3]],
            RA.lst[[4]],
            RA.lst[[5]],
            nrow=5, align = 'v', rel_heights = rep(1,5)/5)

  # Closing the graphical device
  dev.off()

  # Customizing the output
  pdf("./figures/Supp_4_AA_heatmap.pdf",         # File name
      width = 7.67, height = 3.78, # Width and height in inches
      bg = "white")

  # Creating a plot
  plot_grid(AA.lst[[1]],
            AA.lst[[2]],
            AA.lst[[3]],
            AA.lst[[4]],
            AA.lst[[5]],
            nrow=5, align = 'v', rel_heights = rep(1,5)/5)

  # Closing the graphical device
  dev.off()
