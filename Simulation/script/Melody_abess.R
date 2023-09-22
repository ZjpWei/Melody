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
  # arg1: replicate number:
  # from 1 to 100
  s <- as.numeric(args[1])
  
  # arg2: simulation scenario:
  # signature sparsity: "By_taxnum"
  # signature effect direction: "By_pos_prop"
  # signature prevalence: "By_abd_prop"
  # signature effect size: "By_effsz"
  # case/control sequence depth unevenness: "By_seqdepth"
  scenario <- as.numeric(args[2])
  
  # arg3: factor for this scenario: 
  loc <- as.numeric(args[3])

  # parameters
  data.loc <- paste0("./Simulation/", scenario, "/", loc, "/")
  tag <- "Models"
  if(!dir.exists(paste0(data.loc, tag))){
    dir.create(paste0(data.loc, tag))
  }
  ref.1 <- "Coprococcus catus [ref_mOTU_v2_4874]"
  
  # load data
  load(paste0(data.loc, "prepare_data/data.rel.", as.character(s), ".Rdata"))
  data.rel.analysis <- data.rel
  L <- length(data.rel.analysis)

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
  
  # Generate summary statistics
  summary.stat.study <- melody.get.summary(rel.abd = rel.abd,
                                           sample.data = meta,
                                           sample.id = "Sample_ID",
                                           study = "Study",
                                           disease = "Group",
                                           ref = ref.1)
  
  # Meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study, 
                                      tune.path = "sequence",
                                      ouput.best.one = FALSE)
  
  save(Melody.model, file = paste0(data.loc, tag, "/Melody.model.",as.character(s), ".Rdata"))
  