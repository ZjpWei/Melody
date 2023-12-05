# =============================================== #
#           Pooled Meta-Analysis Model            #
# =============================================== #
  
  # Packages ----
  library("miMeta")
  library("ConQuR")
  library("doParallel")
  library("MMUPHin")
  
  # Parameters ----
  cat('Starting model building script\n')
  start.time <- proc.time()[1]
  
  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("The analysis tag needs to be provided! Exiting...\n")
  }
  # arg1: replicate number;
  # from 1 to 100
  s <- as.numeric(args[1])

  # arg2: large/small sample size set;
  # large/samll
  set <- args[2]

  # arg3: simulation scenario:
  # signature sparsity: "Sig_number"
  # signature effect direction: "Sig_effdir"
  # signature effect size: "Sig_effsz"
  # case/control sequence depth unevenness: "Sig_depth"
  scenario <- args[3]

  # arg4: factor for this scenario:
  loc <- args[4]

  # General ----
  data.loc <- paste0("./Simulation/", set, "/", scenario, "/", loc, "/")
  tag <- "Models_original"
  if(!dir.exists(paste0(data.loc, tag))){
    dir.create(paste0(data.loc, tag))
  }
  ref <- "Coprococcus catus [ref_mOTU_v2_4874]"
  
  # Original data ----
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
  meta <- data.frame(Sample_ID = colnames(rel.abd), Pool = as.character(1),
                     Study = Study, Group = Group)
  meta$Study <- factor(meta$Study, levels = c("1", "2", "3", "4", "5"))
  rownames(meta) <- meta$Sample_ID
  
  # Generate summary statistics
  summary.stat.study <- melody.get.summary(rel.abd = rel.abd,
                                           sample.data = meta,
                                           sample.id = "Sample_ID",
                                           study = "Pool",
                                           disease = "Group",
                                           covariates = "Study",
                                           ref = ref)
  
  # Meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study,
                                      tune.path = "sequence",
                                      ouput.best.one = FALSE)
  
  save(Melody.model, file = paste0(data.loc, tag, "/Melody.pool.org.model.",as.character(s), ".Rdata"))

  # MMUPHin batch-corrected ----
  tag <- "Models_MMUPHin_batch_corrected"
  if(!dir.exists(paste0(data.loc, tag))){
    dir.create(paste0(data.loc, tag))
  }
  batch_adj <- MMUPHin::adjust_batch(feature_abd = rel.abd,
                                     batch = "Study",
                                     covariates = "Group",
                                     data = meta)
  
  rel.abd.mmu <- batch_adj$feature_abd_adj
  
  # Generate summary statistics
  summary.stat.study.mmu <- melody.get.summary(rel.abd = rel.abd.mmu,
                                               sample.data = meta,
                                               sample.id = "Sample_ID",
                                               study = "Pool",
                                               disease = "Group",
                                               covariates = "Study",
                                               ref = ref)
  
  # Meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study,
                                      tune.path = "sequence",
                                      ouput.best.one = FALSE)
  
  save(Melody.model, file = paste0(data.loc, tag, "/Melody.pool.mmu.model.",as.character(s), ".Rdata"))
  
  # ConQuR batch-corrected ----
  tag <- "Models_ConQuR_batch_corrected"
  if(!dir.exists(paste0(data.loc, tag))){
    dir.create(paste0(data.loc, tag))
  }
  rel.abd.con <- ConQuR(tax_tab=t(rel.abd), batchid=meta$Study, covariates=meta$Group, 
                        batch_ref="1", num_core = 2)
  
  # Generate summary statistics
  summary.stat.study.bat <- melody.get.summary(rel.abd = t(rel.abd.con),
                                               sample.data = meta,
                                               sample.id = "Sample_ID",
                                               study = "Pool",
                                               disease = "Group",
                                               covariates = "Study",
                                               ref = ref)
  # meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study.bat,
                                      tune.path = "sequence",
                                      ouput.best.one = FALSE)
  
  save(Melody.model, file = paste0(data.loc, tag, "/Melody.pool.con.model.",as.character(s), ".Rdata"))
  
  
