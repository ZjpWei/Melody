# =============================================== #
#                Meta-Analysis Model              #
# =============================================== #

  # Packages ----
  library("miMeta")

  # General ----
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
  # signature sparsity: "Sig_number"
  # signature effect direction: "Sig_effdir"
  # signature effect size: "Sig_effsz"
  # case/control sequence depth unevenness: "Sig_depth"
  scenario <- args[2]

  # arg3: factor for this scenario:
  loc <- args[3]

  # Parameters ----
  data.loc <- paste0("./Simulation/", scenario, "/", loc, "/")
  tag <- "Models_original"
  if(!dir.exists(paste0(data.loc, tag))){
    dir.create(paste0(data.loc, tag))
  }
  ref <- "Coprococcus catus [ref_mOTU_v2_4874]"

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
                                           ref = ref)

  # Meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study,
                                      tune.path = "sequence",
                                      ouput.best.one = FALSE)

  save(Melody.model, file = paste0(data.loc, tag, "/Melody.model.",as.character(s), ".Rdata"))
