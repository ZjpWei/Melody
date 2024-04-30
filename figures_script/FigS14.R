# =============================================== #
#             Supplementary figure 14              #
# =============================================== #

  # Packages ----
  library("miMeta")
  library("ggplot2")
  library("tidyverse")
  library("latex2exp")

  # General ----
  rm(list = ls())
  source("./utility/meta_utility.R")
  load("./CRC/Sensitivity/Melody.model.ref1.Rdata")
  load("./CRC/Sensitivity/summary.stats.ref1.Rdata")

  ## Main script
  delta <- Melody.model.ref1$disease$delta
  best_s <- sum(Melody.model.ref1$disease$coef!=0)
  feature.ID <- NULL
  cov.int.ID <- NULL
  ref <- NULL
  for(d in names(summary.stat.study.ref1)){
    feature.ID <- c(feature.ID, rownames(summary.stat.study.ref1[[d]]$est), summary.stat.study.ref1[[d]]$ref)
    cov.int.ID <- c(cov.int.ID, colnames(summary.stat.study.ref1[[d]]$est))
    ref <- c(ref, summary.stat.study.ref1[[d]]$ref)
  }
  cov.int.ID <- sort(unique(cov.int.ID))
  feature.ID <- sort(unique(feature.ID))
  K <- length(feature.ID)
  
  feature.set <- list()
  for(d in names(summary.stat.study.ref1)){
    taxa.vec.tmp <- rep(FALSE, K)
    names(taxa.vec.tmp) <- feature.ID
    taxa.vec.tmp[rownames(summary.stat.study.ref1[[d]]$est)] <- TRUE
    feature.set[[d]] <- taxa.vec.tmp
  }
  
  summary.stat.study <- list()
  for(d in names(summary.stat.study.ref1)){
    summary.stat.study[[d]] <- list(est = summary.stat.study.ref1[[d]]$est[,1],
                                    cov = diag(summary.stat.study.ref1[[d]]$var[,1]),
                                    n = summary.stat.study.ref1[[d]]$n)
  }
  
  lasso.mat <- Get_lasso_pre(summary.stat.study = summary.stat.study,
                             feature.ID = feature.ID,
                             feature.set = feature.set,
                             study.ID = names(summary.stat.study.ref1),
                             ref = ref,
                             K = K)
  
  ## change s ----
  # Run optimization with shortened support.sizes and best reference taxon estimates
  GIC.s <- NULL
  for(s.lambda in 20:110){
    a <- byabess_AA(summary.stat.study = summary.stat.study, 
                    lasso.mat = lasso.mat,
                    delta = delta, 
                    support.size = s.lambda)
    GIC.s <- c(GIC.s, GIC.cal(result = a, tune.type = "BIC"))
  }
 
  ## change delta 1 ----
  GIC.r1 <- NULL
  ref1.rg <- seq(from = -1, to = 1, length.out = 100)
  tmp.delta <- delta
  for(ref1 in ref1.rg){
    tmp.delta[1] <- ref1
    a <- byabess_AA(summary.stat.study = summary.stat.study, 
                    lasso.mat = lasso.mat,
                    delta = tmp.delta,
                    support.size = best_s)
    GIC.r1 <- c(GIC.r1, GIC.cal(result = a, tune.type = "BIC"))
  }
  
  ## change delta 2 ----
  GIC.r2 <- NULL
  ref2.rg <- seq(from = -1, to = 1, length.out = 100)
  tmp.delta <- delta
  for(ref2 in ref2.rg){
    tmp.delta[2] <- ref2
    a <- byabess_AA(summary.stat.study = summary.stat.study, 
                    lasso.mat = lasso.mat,
                    delta = tmp.delta,
                    support.size = best_s)
    GIC.r2 <- c(GIC.r2, GIC.cal(result = a, tune.type = "BIC"))
  }

  ## change delta 3 ----
  GIC.r3 <- NULL
  ref3.rg <- seq(from = -1, to = 1, length.out = 300)
  tmp.delta <- delta
  for(ref3 in ref3.rg){
    tmp.delta[3] <- ref3
    a <- byabess_AA(summary.stat.study = summary.stat.study, 
                    lasso.mat = lasso.mat,
                    delta = tmp.delta,
                    support.size = best_s)
    GIC.r3 <- c(GIC.r3, GIC.cal(result = a, tune.type = "BIC"))
  }

  ## change delta 4 ----
  GIC.r4 <- NULL
  ref4.rg <- seq(from = -1, to = 1, length.out = 400)
  tmp.delta <- delta
  for(ref4 in ref4.rg){
    tmp.delta[4] <- ref4
    a <- byabess_AA(summary.stat.study = summary.stat.study, 
                    lasso.mat = lasso.mat,
                    delta = tmp.delta,
                    support.size = best_s)
    GIC.r4 <- c(GIC.r4, GIC.cal(result = a, tune.type = "BIC"))
  }
 
  ## change delta 5 ----
  GIC.r5 <- NULL
  ref5.rg <- seq(from = -1, to = 1, length.out = 500)
  tmp.delta <- delta
  for(ref5 in ref5.rg){
    tmp.delta[5] <- ref5
    a <- byabess_AA(summary.stat.study = summary.stat.study, 
                    lasso.mat = lasso.mat,
                    delta = tmp.delta,
                    support.size = best_s)
    GIC.r5 <- c(GIC.r5, GIC.cal(result = a, tune.type = "BIC"))
  }

  # Generate figures ----
  p1 <- data.frame(s = 20:110, BIC = GIC.s) %>% ggplot(aes(x = s, y = BIC)) +
    geom_line() + scale_y_continuous(limits = c(4.2, 4.4), breaks = seq(4.2,4.4, length.out = 3)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  p2 <- data.frame(s = ref1.rg, BIC = GIC.r1) %>% ggplot(aes(x = s, y = BIC)) +
    geom_line() + xlab(TeX('$\\delta_1$')) +
    scale_y_continuous(limits = c(4,7), breaks = seq(4,7, length.out = 7)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  p3 <- data.frame(s = ref2.rg, BIC = GIC.r2) %>% ggplot(aes(x = s, y = BIC)) +
    geom_line() + xlab(TeX('$\\delta_2$')) +
    scale_y_continuous(limits = c(4,7), breaks = seq(4,7, length.out = 7)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))


  p4 <- data.frame(s = ref3.rg, BIC = GIC.r3) %>% ggplot(aes(x = s, y = BIC)) +
    geom_line() + xlab(TeX('$\\delta_3$')) +
    scale_y_continuous(limits = c(4,7), breaks = seq(4,7, length.out = 7)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  p5 <- data.frame(s = ref4.rg, BIC = GIC.r4) %>% ggplot(aes(x = s, y = BIC)) +
    geom_line() + xlab(TeX('$\\delta_4$')) +
    scale_y_continuous(limits = c(4,7), breaks = seq(4,7, length.out = 7)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  p6 <- data.frame(s = ref5.rg, BIC = GIC.r5) %>% ggplot(aes(x = s, y = BIC)) +
    geom_line() + xlab(TeX('$\\delta_5$')) +
    scale_y_continuous(limits = c(4,7), breaks = seq(4,7, length.out = 7)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  ggp1 <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)

  pdf("./figures/FigS14.pdf", width = 13.67, height = 7.29, bg = "white")

  ggp1

  dev.off()
