  ##############################################
#                                                #
#             Supplementary figure 6             #
#                                                #
  ##############################################

  # Packages
  library("miMeta")
  library("ggplot2")
  library("tidyverse")
  library("latex2exp")

  rm(list = ls())
  load("./CRC_meta_analysis/CRC_all_K401/sensitivity/Melody.model.ref1.Rdata")
  best_ref <- Melody.model.ref1$delta
  best_s <- sum(Melody.model.ref1$coef!=0)
  load("./CRC_meta_analysis/CRC_all_K401/summary.stats/summary.stats.ref1.Rdata")
  source("./utility/meta_utility.R")
  summary.stat.study.ref1 <- Get_lasso_pre(Melody = summary.stat.study.ref1)
  #################### change s ####################
  # Run optimzation with  shortened support.sizes, tau and best reference taxon estimates
  GIC.s <- NULL
  for(s.lambda in 20:110){
    a <- byabess_AA(Melody = summary.stat.study.ref1, ref_est = best_ref, support.size = s.lambda)
    GIC.s <- c(GIC.s, GIC.cal(result = a, tune.type = "HBIC"))
  }
  ##############################################
  summary.stat.study <- summary.stat.study.ref1$summary.stat.study
  taxa.set <- summary.stat.study.ref1$taxa.set
  L <- summary.stat.study.ref1$dat.inf$L
  K <- summary.stat.study.ref1$dat.inf$K
  quantile_sm <- NULL
  for(l in 1:L){
    tmp_qt <- quantile(summary.stat.study[[l]]$est, probs = c(0, 1))
    quantile_sm <- rbind(quantile_sm, tmp_qt)
  }
  ################ change delta 1 ##############
  GIC.r1 <- NULL
  ref1.rg <- seq(from = -1, to = 1, length.out = 100)
  tmp.ref <- best_ref
  for(ref1 in ref1.rg){
    tmp.ref[1] <- ref1
    a <- byabess_AA(Melody = summary.stat.study.ref1,
                    ref_est = tmp.ref,
                    support.size = 61)
    GIC.r1 <- c(GIC.r1, GIC.cal(result = a, tune.type = "HBIC"))
  }
  ################ change delta 2 ##############
  GIC.r2 <- NULL
  ref2.rg <- seq(from = -1, to = 1, length.out = 100)
  tmp.ref <- best_ref
  for(ref2 in ref2.rg){
    tmp.ref[2] <- ref2
    a <- byabess_AA(Melody = summary.stat.study.ref1,
                    ref_est = tmp.ref,
                    support.size = 61)
    GIC.r2 <- c(GIC.r2, GIC.cal(result = a, tune.type = "HBIC"))
  }
  ################ change delta 3 ##############
  GIC.r3 <- NULL
  ref3.rg <- seq(from = -1, to = 1, length.out = 300)
  tmp.ref <- best_ref
  for(ref3 in ref3.rg){
    tmp.ref[3] <- ref3
    a <- byabess_AA(Melody = summary.stat.study.ref1,
                    ref_est = tmp.ref,
                    support.size = 61)
    GIC.r3 <- c(GIC.r3, GIC.cal(result = a, tune.type = "HBIC"))
  }
  ################ change delta 4 ##############
  GIC.r4 <- NULL
  ref4.rg <- seq(from = -1, to = 1, length.out = 400)
  tmp.ref <- best_ref
  for(ref4 in ref4.rg){
    tmp.ref[4] <- ref4
    a <- byabess_AA(Melody = summary.stat.study.ref1,
                    ref_est = tmp.ref,
                    support.size = 61)
    GIC.r4 <- c(GIC.r4, GIC.cal(result = a, tune.type = "HBIC"))
  }
  ################ change delta 5 ##############
  GIC.r5 <- NULL
  ref5.rg <- seq(from = -1, to = 1, length.out = 500)
  tmp.ref <- best_ref
  for(ref5 in ref5.rg){
    tmp.ref[5] <- ref5
    a <- byabess_AA(Melody = summary.stat.study.ref1,
                    ref_est = tmp.ref,
                    support.size = 61)
    GIC.r5 <- c(GIC.r5, GIC.cal(result = a, tune.type = "HBIC"))
  }

  # Plotting
  p1 <- data.frame(s = 20:110, HBIC = GIC.s) %>% ggplot(aes(x = s, y = HBIC)) +
    geom_line() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  p2 <- data.frame(s = ref1.rg, HBIC = GIC.r1) %>% ggplot(aes(x = s, y = HBIC)) +
    geom_line() + xlab(TeX('$\\delta_1$')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  p3 <- data.frame(s = ref2.rg, HBIC = GIC.r2) %>% ggplot(aes(x = s, y = HBIC)) +
    geom_line() + xlab(TeX('$\\delta_2$')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))


  p4 <- data.frame(s = ref3.rg, HBIC = GIC.r3) %>% ggplot(aes(x = s, y = HBIC)) +
    geom_line() + xlab(TeX('$\\delta_3$')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  p5 <- data.frame(s = ref4.rg, HBIC = GIC.r4) %>% ggplot(aes(x = s, y = HBIC)) +
    geom_line() + xlab(TeX('$\\delta_4$')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  p6 <- data.frame(s = ref5.rg, HBIC = GIC.r5) %>% ggplot(aes(x = s, y = HBIC)) +
    geom_line() + xlab(TeX('$\\delta_5$')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=16))

  ggp1 <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)

  # Customizing the output
  pdf("./figures/Supp_7.pdf",         # File name
      width = 13.67, height = 7.29, # Width and height in inches
      bg = "white")

  ggp1

  # Closing the graphical device
  dev.off()
