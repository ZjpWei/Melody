  ##############################################
#                                                #
#          Meta-Analysis Model Building          #
#                                                #
  ##############################################

  # Packages
  library('coin')
  library('phyloseq')
  library('pracma')
  library('corpcor')
  library('abess')
  library("miMeta")
  library('MASS')
  library("R6")

  # general
  cat('Starting model building script\n')
  start.time <- proc.time()[1]

  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("The analysis tag needs to be provided! Exiting...\n")
  }
  s <- as.numeric(args[1])

  # Parameters
  data.loc <- "./CRC_meta_analysis/Prediction/Random_split/"
  seed <- 2023
  tag <- "Models_original"
  if(!dir.exists(paste0(data.loc, tag))){
    dir.create(paste0(data.loc, tag))
  }
  ref.1 <- "Coprococcus catus [ref_mOTU_v2_4874]"

  # Load data
  load(paste0(data.loc, "/prepare_data/data.rel.all.Rdata"))
  L <- length(data.rel)
  K <- ncol(data.rel[[1]]$Y)

  # sample case and control
  set.seed(s + seed)
  case.id.lst <- list()
  control.id.lst <- list()
  for(l in 1:L){
    case.id.lst[[l]] <- sample(which(data.rel[[l]]$X == 1))
    control.id.lst[[l]] <- sample(which(data.rel[[l]]$X == 0))
  }

  # Split analysis data, training data and testing data
  data.rel.analysis <- list()
  data.rel.train <- list()
  data.rel.test <- list()
  set.seed(s + seed)
  for(l in 1:L){
    case.id <- case.id.lst[[l]]
    control.id <- control.id.lst[[l]]
    case.rg <- round(c(0.6,0.8,1) * length(case.id))
    control.rg <- round(c(0.6,0.8,1) * length(control.id))
    data.rel.analysis[[l]] <- list(Y = data.rel[[l]]$Y[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]])),],
                                   X = data.rel[[l]]$X[sort(c(case.id[1:case.rg[1]], control.id[1:control.rg[1]]))])

    data.rel.train[[l]] <- list(Y = data.rel[[l]]$Y[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]])),],
                                X = data.rel[[l]]$X[sort(c(case.id[(case.rg[1]+1):case.rg[2]], control.id[(control.rg[1]+1):control.rg[2]]))])

    data.rel.test[[l]] <- list(Y = data.rel[[l]]$Y[sort(c(case.id[(case.rg[2]+1):case.rg[3]], control.id[(control.rg[2]+1):control.rg[3]])),],
                               X = data.rel[[l]]$X[sort(c(case.id[(case.rg[2]+1):case.rg[3]], control.id[(control.rg[2]+1):control.rg[3]]))])
  }

  Study <- NULL
  Group <- NULL
  rel.abd <- NULL
  for(l in 1:length(data.rel.analysis)){
    Study <- c(Study, rep(as.character(l), nrow(data.rel.analysis[[l]]$Y)))
    Group <- c(Group, as.character(data.rel.analysis[[l]]$X))
    rel.abd <- rbind(rel.abd, data.rel.analysis[[l]]$Y)
  }
  rel.abd <- t(rel.abd)
  meta <- data.frame(Sample_ID = colnames(rel.abd), Study = Study, Group = Group)
  rownames(meta) <- meta$Sample_ID

  # generate summary statistics
  summary.stat.study <- melody.get.summary(rel.abd = rel.abd,
                                           sample.data = meta,
                                           sample.id = "Sample_ID",
                                           study = "Study",
                                           disease = "Group",
                                           ref = ref.1,
                                           verbose = TRUE)

  # meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study, tune.path = "sequence", tune.size.sequence = c(1:100))

  save(Melody.model, file = paste0(data.loc, tag, "/Melody.model.",as.character(s),".Rdata"))
