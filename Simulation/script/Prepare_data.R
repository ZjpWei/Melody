# =============================================== #
#            Data preparation script              #
# =============================================== #

  # Packages ----
  library('phyloseq')
  library('coin')
  library("purrr")

  cat('Starting data cleaning script\n')
  start.time <- proc.time()[1]
  # Parameters ----
  # Set active taxa, the number should be Ka
  signal.idx <- NA
  if(!is.na(signal.idx)){
    stopifnot(is.numeric(signal.idx) & is.vector(signal.idx) & length(signal.idx  == Ka))
  }

  # Global random seed
  seed <- 2023
  # General ----
  args = commandArgs(trailingOnly=TRUE)
  if(length(args)==0){
    stop("The analysis tag needs to be provided! Exiting...\n")
  }
  print(args)
  # replicate id in simulation
  s <- as.numeric(args[1])
  # signature sparsity {20, 30, 40(default), 50, 60, 70, 80}
  Ka <- as.numeric(args[2])
  # signature effect direction {0.5, 0.6, 0.7(default), 0.8, 0.9, 1}
  pos.pt <- as.numeric(args[3])
  # signature effect size {Lager: 1, 1.5, 2(default), 2.5, 3
  #                        Small: 4, 4.5, 5(default), 5.5, 6}
  effect.sz <- as.numeric(args[4])
  # case/control sequence depth unevenness {0(default), 0.25, 0.5, 0.75, 1}
  mu <- as.numeric(args[5])
  # Large/Small sample set {small, large}
  set <- as.numeric(args[6])
  
  # set document location respecting the input parameter above
  # signature sparsity: "Sig_number"
  # signature effect direction: "Sig_effdir"
  # signature effect size: "Sig_effsz"
  # case/control sequence depth unevenness: "Sig_depth"

  # Sample size vector in absolute data simulation
  if(set == "large"){
    n.sample <- c(100, 120, 140, 160, 180)
    loc.vec <- c(Ka == 40, pos.pt == 0.7, effect.sz == 2, mu == 0)
  }else if(set == "small"){
    n.sample <- c(20, 30, 40, 50, 60)
    loc.vec <- c(Ka == 40, pos.pt == 0.7, effect.sz == 5, mu == 0)
  }
  if(sum(loc.vec) < 3){
    stop("Please only change one factor at the same time.\n")
  }else if(sum(loc.vec) == 4){
    data.loc <- c(paste0("./Simulation/", set, "/Sig_number/d", as.character(Ka), "/"),
                  paste0("./Simulation/", set, "/Sig_effdir/p", as.character(pos.pt * 100), "/"),
                  paste0("./Simulation/", set, "/Sig_effsz/sig", as.character(effect.sz), "/"),
                  paste0("./Simulation/", set, "/Sig_depth/seq", as.character(mu), "/"))
  }else{
    if(which(!loc.vec) == 1){
      data.loc <-  paste0("./Simulation/", set, "/Sig_number/d", as.character(Ka), "/")
    }else if(which(!loc.vec) == 2){
      data.loc <-  paste0("./Simulation/", set, "/Sig_effdir/p", as.character(pos.pt * 100), "/")
    }else if(which(!loc.vec) == 3){
      data.loc <-  paste0("./Simulation/", set, "/Sig_effsz/sig", as.character(effect.sz), "/")
    }else if(which(!loc.vec) == 4){
      data.loc <-  paste0("./Simulation/", set, "/Sig_depth/seq", as.character(mu), "/")
    }
  }
  for(loc in data.loc){
    if(!dir.exists(loc)){
      dir.create(loc)
      dir.create(paste0(loc, "prepare_data"))
      dir.create(paste0(loc, "signal"))
    }else{
      if(!dir.exists(paste0(loc, "prepare_data"))){
        dir.create(paste0(loc, "prepare_data"))
      }
      if(!dir.exists(paste0(loc, "signal"))){
        dir.create(paste0(loc, "signal"))
      }
    }
  }
  
  # load utility
  source("./utility/GDM_utility.R")

  # load count data
  load("./Simulation/data/count.Rdata")

  # load meta data
  load("./Simulation/data/meta.Rdata")

  # data processing
  abd.pt <- 0.5
  meta <- as.data.frame(meta)
  study <- c("AT-CRC", "CN-CRC", "DE-CRC", "FR-CRC", "US-CRC")
  rownames(meta) <- meta$Sample_ID
  meta$Group <- factor(meta$Group, level = c("CTR", "CRC"))
  meta$Study <- factor(meta$Study, levels = study)
  sample.id.kp <- names(which(rowSums(count) >= 2000))
  meta <- meta[sample.id.kp,]
  count <- count[sample.id.kp,]
  len.filter <- colMeans(count != 0) >= 0.2
  count <- count[,len.filter]
  
  # split data to corresponding study
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
  
  # set random seed
  set.seed(s + seed)
  
  # simulate null data
  data.null <- list()
  for(l in 1:L){
    load(paste0("./Simulation/GDM_fit/GDM.", as.character(l), ".Rdata"))
    X <- matrix(1, nrow = n.sample[l], ncol = 1)
    dmData <- r.GDM(X = X, mod = mod.gdm)
    Y.tmp <- Y.study[[l]]
    idx <- order(colMeans(Y.tmp/rowSums(Y.tmp)), decreasing = TRUE)
    idx.rev <- rank(-colMeans(Y.tmp/rowSums(Y.tmp)))
    stopifnot(sum(colMeans(Y.tmp/rowSums(Y.tmp)) != 0) == ncol(dmData))
    dmData <- cbind(dmData, matrix(0, nrow = n.sample[l], ncol = K - ncol(dmData)))[,idx.rev]
    colnames(dmData) <- colnames(Y.study[[l]])
    data.null[[l]] <- list(Y = dmData, X = rep(0, nrow(dmData)))
  }
  
  # decide the active taxa
  if(is.na(signal.idx)){
    # random signal
    ave.prop <- colMeans(count / rowSums(count))
    # decide most abundant taxa and less abundant taxa
    signal.most.abund <- which(ave.prop >= 1e-3)
    signal.less.abund <- which(ave.prop  < 1e-3)
    signal.m.abd <- sample(signal.most.abund)
    signal.l.abd <- sample(signal.less.abund)
    tmp.signal.idx <- c(signal.m.abd[1:round(abd.pt * Ka)], signal.l.abd[1:(Ka - round(abd.pt * Ka))])
    # decide signal location and direction
    signal.sig <- sample(c(rep(TRUE, round(pos.pt * Ka)), rep(FALSE, Ka - round(pos.pt*Ka))))
    signal.idx = sort(tmp.signal.idx)
  }
  
  # add signal
  signal.add <- c()
  tmp.add <- runif(n = Ka, min = 1, max = effect.sz + 1)
  for(l in 1:L){
    signal.add <- rbind(signal.add, tmp.add)
  }
  
  # Add signal to null data
  # simulate the absolute abundant data and compute corresponded relative aboudant data
  data.rel <- list()
  for(l in 1:L){
    org.data <- Y.study[[l]]
    c.name <- colnames(org.data)
    n = n.sample[l]
    # generate absolute abundant data 
    # shuffle subjects (we always assign the first half subjects to cases)
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
    
    # generate relative abundant data 
    Simulate.depth.1 <- rep(0, n)
    Simulate.depth.1[case.idx.1] <- sample(x = rowSums(Y.study[[l]]), size = length(case.idx.1), replace = TRUE) * (mu + 1) 
    Simulate.depth.1[-case.idx.1]<- sample(x = rowSums(Y.study[[l]]), size = n - length(case.idx.1), replace = TRUE)
    Simulate.count.rel.1 = NULL
    for(ll in 1:nrow(Simulate.count.1)){
      if(Simulate.otc.1[ll] == 1){
        Simulate.count.rel.1 <- cbind(Simulate.count.rel.1, rmultinom(1, Simulate.depth.1[ll], Prob.abs.1[ll,]))
      }else{
        Simulate.count.rel.1 <- cbind(Simulate.count.rel.1, rmultinom(1, Simulate.depth.1[ll], Prob.abs.1[ll,]))
      }
    }
    Simulate.count.rel.1 <- t(Simulate.count.rel.1)
    # generate taxa and sample names
    data.rel[[l]] <- list(Y=Simulate.count.rel.1, X=Simulate.otc.1)
    sample.nm <- paste0("Cohort.",as.character(l),".s",as.character(1:n))
    colnames(data.rel[[l]]$Y) <- c.name
    rownames(data.rel[[l]]$Y) <- sample.nm
  }
  
  # data summary ----
  # cases and controls
  signal.idx <- colnames(data.rel[[1]]$Y)[signal.idx]
  for(loc in data.loc){
    save(data.rel, file = paste0(loc, "prepare_data/data.rel.", as.character(s),".Rdata"))
    
    save(signal.idx, file = paste0(loc, "signal/signal.", as.character(s),".Rdata"))
  }
  
