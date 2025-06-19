# =============================================== #
#             Supplementary figure 6              #
# =============================================== #

  # Packages ----
  library("pheatmap")
  library("miMeta")

  # General ----
  rm(list = ls())
  load("./CRC/Processed_data/data.org.K401.Rdata")
  L <- length(data.rel)

  # Train Melody model with sequence:
  L <- length(data.rel)
  rel.abd <- list()
  refs <- c()
  covariate.interest <- list()
  for(l in 1:length(data.rel)){
    rel.abd[[paste0("S",l)]] <- data.rel[[l]]$Y
    covariate.interest[[paste0("S",l)]] <- data.frame("disease" = data.rel[[l]]$X)
    refs <- c(refs, "Coprococcus catus [ref_mOTU_v2_4874]")
  }
  names(refs) <- names(covariate.interest)

  null.obj <- melody.null.model(rel.abd = rel.abd, ref = refs)

  summary.stat.study <- melody.get.summary(null.obj = null.obj,
                                           covariate.interest = covariate.interest)

  # Meta-analysis
  Melody.model <- melody.meta.summary(summary.stats = summary.stat.study,
                                      tune.path = "sequence",
                                      tune.size.sequence = 1:60,
                                      output.best.one = FALSE,
                                      verbose = TRUE)

  load("./Data/CRC_data/data/tax.Rdata")
  load("./CRC/Models_org/BW.prop.model.Rdata")
  load("./CRC/Models_org/ANCOMBC2.model.Rdata")
  load("./CRC/Models_org/Aldex2.model.Rdata")
  load("./CRC/Models_org/MMUPHin.model.Rdata")
  load("./CRC/Models_org/CLR.lasso.model.Rdata")

  # Generate figures ----
  for(tops in c(20, 30, 40)){
    index <- list()
    taxa.id <- gsub("[][]", "",unlist(regmatches(names(which(Melody.model$disease$coef[,tops]!=0)),
                                                 gregexpr("\\[.*?\\]",names(which(Melody.model$disease$coef[,tops]!=0))))))
    len.tax <- length(taxa.id)
    index[[1]] <- taxa.id
    ###
    aldex2.top <- rownames(Aldex2.model$glmfit)[rank(Aldex2.model$glmfit[,c("labels1:pval")]) <= len.tax]
    taxa.id <- gsub("[][]", "",unlist(regmatches(aldex2.top, gregexpr("\\[.*?\\]",aldex2.top))))
    index[[2]] <- taxa.id
    ###
    ancom.top <- ANCOMBC2.model$res$taxon[rank(ANCOMBC2.model$res$p_labels1) <= len.tax]
    taxa.id <- gsub("[][]", "",unlist(regmatches(ancom.top, gregexpr("\\[.*?\\]",ancom.top))))
    index[[3]] <- taxa.id
    ###
    mmuphin.top <- MMUPHin.model$meta_fits$feature[rank(MMUPHin.model$meta_fits$pval) <= len.tax]
    taxa.id <- gsub("[][]", "",unlist(regmatches(mmuphin.top, gregexpr("\\[.*?\\]",mmuphin.top))))
    index[[4]] <- taxa.id
    ###
    BW.top <- rownames(BW.prop.model)[rank(BW.prop.model$p.val) <= len.tax]
    taxa.id <- gsub("[][]", "",unlist(regmatches(BW.top, gregexpr("\\[.*?\\]",BW.top))))
    index[[5]] <- taxa.id
    ###
    len.tmp <- c()
    for(ll in 1:ncol(clrlasso$glmnet.fit$beta)){
      len.tmp <- c(len.tmp, sum(clrlasso$glmnet.fit$beta[,ll]!=0))
    }
    taxa.id <- names(which(rank(-abs(clrlasso$glmnet.fit$beta[,which(len.tmp >= len.tax)[1]])) <= len.tax))
    taxa.id <- gsub("[][]", "",unlist(regmatches(taxa.id, gregexpr("\\[.*?\\]",taxa.id))))
    index[[6]] <- taxa.id

    # prepare plotting ----
    index <- index[c(6,3,2,5,1,4)]
    Heatmap <- matrix(0,6,6)
    for(l in 1:length(index)){
      for(ll in 1:length(index)) {
        Heatmap[l, ll] <- length(intersect(index[[l]], index[[ll]]))/
          length(unique(c(index[[l]], index[[ll]])))
      }
    }
    rownames(Heatmap) <- c("CLR-LASSO", "ANCOM-BC2", "ALDEx2", "BW", "Melody", "MMUPHin")

    pdf(paste0("./figures/FigS6_top", tops,".pdf"), width = 6, height = 5, bg = "white")

    pheatmap::pheatmap(Heatmap, cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(20:100)/100,
                       color = (colorRampPalette(c("#2A788EFF", "#7AD151FF", "#FDE725FF"))(100))[20:100],
                       display_numbers = TRUE,fontsize_number = 20, fontsize = 15, legend = FALSE)

    dev.off()
  }
