  ##############################################
#                                                #
#          Figure 2 (a) Melody heatmap           #
#                                                #
  ##############################################

  # Packages
  library("RColorBrewer")
  library("stringr")
  library("tidyverse")
  library("ggtext")
  library("cowplot")
  library("grid")
  library("gridExtra")
  library("ggplot2")

  rm(list = ls())
  load("./CRC_meta_analysis/data/tax.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/summary.stats/summary.stats.ref1.Rdata")
  load("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.ref1.Rdata")
  summary.stat.study <- summary.stat.study.ref1
  Melody.model <- Melody.model.ref1
  len.melody <- sum(Melody.model$coef!=0)

  # Main
  load("./CRC_meta_analysis/CRC_all_K401/prepare_data/data.rel.all.Rdata")
  L <- length(data.rel)
  species.name <- colnames(data.rel[[1]]$Y)
  tab.AA <- NULL
  for(l in 1:L){
    tmp.AA <- rep(NA, ncol(data.rel[[1]]$Y))
    names(tmp.AA) <- colnames(data.rel[[1]]$Y)
    tmp.AA[names(summary.stat.study$summary.stat.study[[l]]$est)] <- summary.stat.study$summary.stat.study[[l]]$est
    tmp.AA[names(Melody.model$delta)[l]] <- 0
    tmp.AA <- tmp.AA + Melody.model$delta[l]
    tab.AA <- rbind(tab.AA, tmp.AA)
  }
  est.vals <- data.frame(all = Melody.model$coef)
  rownames(tab.AA) <- paste0("CRC", as.character(1:L))
  colnames(tab.AA) <- species.name
  tab.AA <- as.data.frame(t(tab.AA))
  ref.studies <- paste0("CRC", as.character(1:L))
  species.heatmap <- rownames(est.vals)[which(est.vals$all != 0)]
  est.vals.signed <- est.vals[species.heatmap,"all", drop=FALSE]
  species.heatmap.orderd <- c(rownames(est.vals.signed[order(est.vals.signed[,"all"],decreasing = T),,drop=FALSE]))

  # Prepare plotting
  est.vals.plot <- est.vals[species.heatmap.orderd, "all"]
  tab.AA.plot <- tab.AA[species.heatmap.orderd,]
  mx <- 2.5
  brs = seq(-mx, mx, by=1e-3)
  num.col.steps = length(brs)-1
  n = floor(0.5*num.col.steps)
  col.hm = c(rev(colorRampPalette(brewer.pal(9, 'Blues'))(n)),
             rep('#FFFFFF', num.col.steps-2*n),
             colorRampPalette(brewer.pal(9, 'Reds'))(n))

  tab.AA.plot[tab.AA.plot >= mx] <- mx
  tab.AA.plot[tab.AA.plot <= -mx] <- -mx

  # Rename the taxa name in plot
  rownames(tax) <- tax$OTUID
  renames <- sub("^\\S+\\s+",
                 '', tax[gsub("[][]", "",unlist(regmatches(species.heatmap.orderd,
                 gregexpr("\\[.*?\\]",species.heatmap.orderd)))),"mOTU"])
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
  renames <- paste0(renames, " (",sub("]",")", sub(".*v2_", "", species.heatmap.orderd),")"))
  renames <- factor(renames, levels=renames)

  # Heatmap plot
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
      geom_tile() + theme_minimal() + ylab("Studies") +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 15),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_rect(fill=NULL, colour='black'),
            plot.margin = unit(c(0, 0, 0, 0), 'cm'),
            axis.text.y = element_text(size = 13),
            axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1, size = 12),
            legend.title = element_text(hjust = 0.5),
            legend.position = "bottom") +
      scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx),
                           breaks = c(-mx,-1,0,1,mx),
                           labels = c(expression(""<="-2.5"), as.character(-1:1), expression("">="2.5"))) +
      guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5, title = "AA coefficient",
                                    title.position = "top"))
  }
  g1 <- plot.single.study.heatmap.ref(rev(ref.studies))

  # Scatter plot
  g2 <- tibble(species=renames,
               est.vals= est.vals.plot,
               variable = "AA beta_hat est",
               colour=est.vals != 0) %>%
    ggplot() + geom_point(aes(x= renames, y= est.vals)) +
    theme_classic() +
    ylab('Meta') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_blank(),
          legend.position = "right") +
    scale_y_continuous(limits=c(min(est.vals.plot) - 0.5, max(est.vals.plot) + 0.5), expand = c(0, 0)) +
    scale_x_discrete(position='bottom') + #facet_wrap(~variable) +
    scale_fill_manual(values=c('lightgrey', 'darkgrey'), guide="none") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed")


   # y.grob <- textGrob(paste0("                                                        ",
   #                          "AA association coefficient"),  gp=gpar(fontsize=15), rot=90)

  pp1 <- plot_grid(g2, g1, nrow=2, align = 'v', rel_heights = c(.2,.8))

  # Customizing the output
  pdf("./figures/figure2_a.pdf",         # File name
      width = 12, height = 5.92, # Width and height in inches
      bg = "white")

  # Creating a plot
  pp1

  # Closing the graphical device
  dev.off()
