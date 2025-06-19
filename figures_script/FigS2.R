# =============================================== #
#             Supplementary figure 2              #
# =============================================== #

  # Packages ----
  library("phyloseq")
  library("ggplot2")
  library("cowplot")
  library("dplyr")

  # General ----
  rm(list = ls())

  ## Original ----
  load("./Data/CRC_data/data/meta.Rdata")
  load("./CRC/Processed_data/data.org.K849.Rdata")
  title <- "CRC Real"
  studymeta <- unique(meta$Study)
  feature.table <- NULL
  Group <- NULL
  Study <- NULL
  K <- 849
  for(l in 1:length(data.rel)){
    feature.table <- rbind(feature.table, data.rel[[l]]$Y/rowSums(data.rel[[l]]$Y) )
    Group <- c(Group, data.rel[[l]]$X)
    Study <- c(Study, rep(paste0("CRC",as.character(l)), length(data.rel[[l]]$X)))
  }
  names(Group) <- rownames(feature.table)
  meta.data = data.frame(labels = Group)
  rownames(feature.table) <- rownames(meta.data)

  feature.table = t(feature.table)
  OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
  META = sample_data(meta.data)
  PHYSEQ = phyloseq(OTU, META)

  # Calculate ordination
  iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "bray")

  all(rownames(iMDS$vectors) == rownames(meta.data))
  PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"],
                     PCo2 = iMDS$vectors[,"Axis.2"],
                     Study = Study,
                     Group = as.character(Group), variable = title)


  ##############################################################################
  ## NetMoss plot
  axis.1.title <- "PCo1"
  axis.2.title <- "PCo2"

  # main plot
  g.main <- PCoA %>%
    ggplot(aes(x=PCo1, y=PCo2, col=Study, shape = Group)) +
    geom_point() +
    scale_shape_manual(breaks = c(0,1), values = c(1, 19)) +
    scale_x_continuous(position='top') +
    xlab(axis.1.title) + ylab(axis.2.title) +
    theme(panel.background = element_rect(fill = 'white'),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.text = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid = element_blank(),
          legend.position = "none")

  #g.main
  # study boxplot axis 1
  g.s.1 <- PCoA %>% mutate(Study=Study) %>%
    ggplot(aes(y=PCo1, x=Study, fill=Study)) +
    geom_boxplot(size=0.1) +
    theme(axis.ticks = element_blank(),
          panel.background = element_rect(fill='white', color = 'black'),
          axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          panel.grid = element_blank(),
          legend.position = "none") +
    coord_flip()

  g.g.1 <- PCoA %>% mutate(Group=Group) %>%
    ggplot(aes(y=PCo1, x=Group, fill=Group)) +
    geom_boxplot(size=0.1) +
    scale_fill_manual(breaks = c(0,1), values = c("grey", "orange"))+
    theme(axis.ticks.y = element_blank(),
          panel.background = element_rect(fill='white', color = 'black'),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          panel.grid = element_blank(),
          legend.position = "none") +
    coord_flip()

  # study boxplot axis 2
  g.s.2 <- PCoA %>% mutate(Study=Study) %>%
    ggplot(aes(y=PCo2, x=Study, fill=Study)) +
    geom_boxplot(size=0.1) +
    scale_y_continuous(position='right') +
    scale_x_discrete(position='top') +
    theme(axis.ticks=element_blank(),
          panel.background = element_rect(fill='white', color = 'black'),
          axis.text = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")

  g.g.2 <- PCoA %>% mutate(Group=Group) %>%
    ggplot(aes(y=PCo2, x=Group, fill=Group)) +
    geom_boxplot(size=0.1) + scale_fill_manual(breaks = c(0,1), values = c("grey", "orange"))+
    scale_y_continuous(position='right') +
    scale_x_discrete(position='top') +
    theme(axis.ticks.x=element_blank(),
          panel.background = element_rect(fill='white', color = 'black'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")


  p.values <- NULL
  w.test <- kruskal.test(PCo1 ~ Study, data = PCoA)
  p.values <- c(p.values, w.test$p.value)
  w.test <- kruskal.test(PCo2 ~ Study, data = PCoA)
  p.values <- c(p.values, w.test$p.value)
  w.test <- kruskal.test(PCo1 ~ Group, data = PCoA)
  p.values <- c(p.values, w.test$p.value)
  w.test <- kruskal.test(PCo2 ~ Group, data = PCoA)
  p.values <- c(p.values, w.test$p.value)
  names(p.values) <- c("Study PCo1", "Study PCo2", "Group PCo1", "Group PCo2")

  # ##############################################################################
  ### Plot  fig 1b
  p.rel <- plot_grid(g.main, g.s.2, g.g.2,
                     g.s.1, NULL, NULL,
                     g.g.1, NULL, NULL, nrow=3, rel_widths = c(0.8, 0.2, 0.13), rel_heights = c(0.8, 0.2, 0.13))

  ### Generate figures ----
  pdf("./figures/FigS2_rel.pdf", width = 9.00, height = 8.30, bg = "white")

  p.rel

  dev.off()

  ## Simulation data ----
  rm(list = ls())
  set.seed(2025)
  n.sample <- c(100, 120, 140, 160, 180)
  Ka <- 40
  pos.pt <- 0.6
  effect.sz <- 2
  mu <-  0
  abd.pt <- 0.5

  ## Data processing, remove samples with sequence depth < 2000, apply 0.2 prevalence filter.
  load("./Data/CRC_data/data/count.Rdata")
  load("./Data/CRC_data/data/meta.Rdata")

  meta <- as.data.frame(meta)
  study <- c("AT-CRC", "CN-CRC", "DE-CRC", "FR-CRC", "US-CRC")
  rownames(meta) <- meta$Sample_ID
  meta$Group <- factor(meta$Group, level = c("CTR", "CRC"))
  meta$Study <- factor(meta$Study, levels = study)
  sample.id.kp <- names(which(rowSums(count) >= 2000))
  meta <- meta[sample.id.kp,]
  count <- count[sample.id.kp,]
  pre.filter <- colMeans(count != 0) >= 0.2
  count <- count[,pre.filter]

  L <- length(study)
  K <- ncol(count)
  Y.study <- list()
  X.study <- list()
  for(l in 1:L){
    id.set <- meta$Study == study[l]
    condition <- meta$Group[id.set]
    Y.study[[l]] <- count[id.set,]
    X.study[[l]] <- c(condition == 'CRC') + 0
  }

  ## Generate signal and signature size
  source("./utility/GDM_utility.R")
  set.seed(1234)
  ## Randomize uneveness
  uneveness.ind <- rbernoulli(n = L, p = 0.5)
  ## Randomize signal
  ave.prop <- colMeans(count / rowSums(count))
  ## Detect most abundant taxa and less abundant taxa
  signal.most.abund <- which(ave.prop >= 1e-3)
  signal.less.abund <- which(ave.prop  < 1e-3)
  signal.m.abd <- sample(signal.most.abund)
  signal.l.abd <- sample(signal.less.abund)
  tmp.signal.idx <- c(signal.m.abd[1:round(abd.pt * Ka)], signal.l.abd[1:(Ka - round(abd.pt * Ka))])
  ## Decide signal location and direction
  signal.sig <- sample(c(rep(TRUE, round(pos.pt * Ka)), rep(FALSE, Ka - round(pos.pt*Ka))))
  signal.idx = sort(tmp.signal.idx)
  ## Add signal
  signal.add <- c()
  tmp.add <- runif(n = Ka, min = 1, max = effect.sz + 1)
  for(l in 1:L){
    signal.add <- rbind(signal.add, tmp.add)
  }

  ## Simulate null AA data
  data.null <- list()
  for(l in 1:L){
    load(paste0("./Simulation/GDM_fit/GDM.", as.character(l), ".Rdata"))
    X <- matrix(1, nrow = n.sample[l], ncol = 1)
    dmData <- r.GDM(X = X, mod = mod.gdm)
    Y.tmp <- Y.study[[l]]
    idx <- order(colMeans(Y.tmp/rowSums(Y.tmp)), decreasing = TRUE)
    idx.rev <- rank(-colMeans(Y.tmp/rowSums(Y.tmp)))
    ## Stop if dimension doesn't match
    stopifnot(sum(colMeans(Y.tmp/rowSums(Y.tmp)) != 0) == ncol(dmData))
    dmData <- cbind(dmData, matrix(0, nrow = n.sample[l], ncol = K - ncol(dmData)))[,idx.rev]
    colnames(dmData) <- colnames(Y.study[[l]])
    data.null[[l]] <- list(Y = dmData, X = rep(0, nrow(dmData)))
  }

  ## Simulate relative abundant data
  data.rel <- list()
  for(l in 1:L){
    org.data <- Y.study[[l]]
    c.name <- colnames(org.data)
    n = n.sample[l]

    ## Shuffle subjects (we always regard the first half subjects as cases)
    case.idx.1 = 1:round(n/2)
    Simulate.count.1 = data.null[[l]]$Y
    Simulate.otc.1 = data.null[[l]]$X
    for(ll in 1:Ka){
      if(signal.sig[ll]){
        Simulate.count.1[case.idx.1, signal.idx[ll]] <- Simulate.count.1[case.idx.1, signal.idx[ll]] * signal.add[l,ll]
      }else{
        Simulate.count.1[-case.idx.1, signal.idx[ll]]<- Simulate.count.1[-case.idx.1,signal.idx[ll]] * signal.add[l,ll]
      }
    }
    Simulate.otc.1[case.idx.1] = 1
    Prob.abs.1 = Simulate.count.1/rowSums(Simulate.count.1)

    ## Simulate relative abundant data
    Simulate.depth.1 <- rep(0, n)
    if(uneveness.ind[l]){
      Simulate.depth.1[case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = length(case.idx.1), replace = TRUE) * (mu + 1)
      Simulate.depth.1[-case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = n - length(case.idx.1), replace = TRUE)
    }else{
      Simulate.depth.1[case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = length(case.idx.1), replace = TRUE)
      Simulate.depth.1[-case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = n - length(case.idx.1), replace = TRUE) * (mu + 1)
    }
    Simulate.count.rel.1 = NULL
    for(ll in 1:nrow(Simulate.count.1)){
      if(Simulate.otc.1[ll] == 1){
        Simulate.count.rel.1 <- cbind(Simulate.count.rel.1, rmultinom(1, Simulate.depth.1[ll], Prob.abs.1[ll,]))
      }else{
        Simulate.count.rel.1 <- cbind(Simulate.count.rel.1, rmultinom(1, Simulate.depth.1[ll], Prob.abs.1[ll,]))
      }
    }
    Simulate.count.rel.1 <- t(Simulate.count.rel.1)

    ## Assign taxa and sample names
    data.rel[[l]] <- list(Y=Simulate.count.rel.1, X=Simulate.otc.1)
    sample.nm <- paste0("Cohort.",as.character(l),".s",as.character(1:n))
    colnames(data.rel[[l]]$Y) <- c.name
    rownames(data.rel[[l]]$Y) <- sample.nm
  }

  feature.table <- NULL
  Group <- NULL
  Study <- NULL
  for(l in 1:length(data.rel)){
    feature.table <- rbind(feature.table, data.rel[[l]]$Y/rowSums(data.rel[[l]]$Y))
    Group <- c(Group, data.rel[[l]]$X)
    Study <- c(Study, rep(paste0("CRC",as.character(l)), length(data.rel[[l]]$X)))
  }
  rownames(feature.table) <- paste0("sa", as.character(1:length(Group)))
  names(Group) <- rownames(feature.table)
  meta.data = data.frame(labels = Group)
  rownames(feature.table) <- rownames(meta.data)

  feature.table = t(feature.table)
  OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
  META = sample_data(meta.data)
  PHYSEQ = phyloseq(OTU, META)

  # Calculate ordination
  iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "bray")

  PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"],
                     PCo2 = iMDS$vectors[,"Axis.2"],
                     Study = Study,
                     Group = as.character(Group))

  ##############################################################################
  ## NetMoss plot
  axis.1.title <- "PCo1"
  axis.2.title <- "PCo2"

  # main plot
  g.main <- PCoA %>%
    ggplot(aes(x=PCo1, y=PCo2, col=Study, shape = Group)) +
    geom_point() +
    scale_shape_manual(breaks = c(0,1), values = c(1, 19)) +
    scale_x_continuous(position='top') +
    xlab(axis.1.title) + ylab(axis.2.title) +
    theme(panel.background = element_rect(fill = 'white'),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.text = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid = element_blank(),
          legend.position = "none")

  #g.main
  # study boxplot axis 1
  g.s.1 <- PCoA %>% mutate(Study=Study) %>%
    ggplot(aes(y=PCo1, x=Study, fill=Study)) +
    geom_boxplot(size=0.1) +
    theme(axis.ticks = element_blank(),
          panel.background = element_rect(fill='white', color = 'black'),
          axis.text = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          panel.grid = element_blank(),
          legend.position = "none") +
    coord_flip()

  g.g.1 <- PCoA %>% mutate(Group=Group) %>%
    ggplot(aes(y=PCo1, x=Group, fill=Group)) +
    geom_boxplot(size=0.1) +
    scale_fill_manual(breaks = c(0,1), values = c("grey", "orange"))+
    theme(axis.ticks.y = element_blank(),
          panel.background = element_rect(fill='white', color = 'black'),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          panel.grid = element_blank(),
          legend.position = "none") +
    coord_flip()

  # study boxplot axis 2
  g.s.2 <- PCoA %>% mutate(Study=Study) %>%
    ggplot(aes(y=PCo2, x=Study, fill=Study)) +
    geom_boxplot(size=0.1) +
    scale_y_continuous(position='right') +
    scale_x_discrete(position='top') +
    theme(axis.ticks=element_blank(),
          panel.background = element_rect(fill='white', color = 'black'),
          axis.text = element_blank(),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")

  g.g.2 <- PCoA %>% mutate(Group=Group) %>%
    ggplot(aes(y=PCo2, x=Group, fill=Group)) +
    geom_boxplot(size=0.1) + scale_fill_manual(breaks = c(0,1), values = c("grey", "orange"))+
    scale_y_continuous(position='right') +
    scale_x_discrete(position='top') +
    theme(axis.ticks.x=element_blank(),
          panel.background = element_rect(fill='white', color = 'black'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")


  p.values <- NULL
  w.test <- kruskal.test(PCo1 ~ Study, data = PCoA)
  p.values <- c(p.values, w.test$p.value)
  w.test <- kruskal.test(PCo2 ~ Study, data = PCoA)
  p.values <- c(p.values, w.test$p.value)
  w.test <- kruskal.test(PCo1 ~ Group, data = PCoA)
  p.values <- c(p.values, w.test$p.value)
  w.test <- kruskal.test(PCo2 ~ Group, data = PCoA)
  p.values <- c(p.values, w.test$p.value)
  names(p.values) <- c("Study PCo1", "Study PCo2", "Group PCo1", "Group PCo2")

  # ##############################################################################
  ### Plot  fig 1b
  #pdf("fig1b.pdf",width = 5,height = 5)

  p.sim <- plot_grid(g.main, g.s.2, g.g.2,
                     g.s.1, NULL, NULL,
                     g.g.1, NULL, NULL, nrow=3, rel_widths = c(0.8, 0.2, 0.13), rel_heights = c(0.8, 0.2, 0.13))

  ### Generate figures ----
  pdf("./figures/FigS2_sim.pdf", width = 9.00, height = 8.30, bg = "white")

  p.sim

  dev.off()
