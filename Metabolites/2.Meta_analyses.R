# =============================================== #
#            Metabolites meta-analysis            #
# =============================================== #  

  rm(list = ls())
  # Package
  library("miMeta")
  library("tidyverse")
  
  # Load data
  load("./Metabolites/Processed_data/processed_data.Rdata")
  cluster <- lapply(metadata, function(d){
    d <- d %>% dplyr::transmute(Sample, Subject) %>% data.frame() 
    rownames(d) <- NULL
    d <- d %>% tibble::column_to_rownames("Sample")
    return(d)
  })
  Metabolite.ids <- unique(common.pairs$Compound)
  Genera.id <- unique(common.pairs$Taxon)
  
  covariates.disease <- lapply(metadata[datasets], function(d){
    case.names <- c("Healthy", "Control", "Normal", "H", "0", "nonIBD", "NA")
    return(d %>% tibble::rownames_to_column(var = "IDS") %>%
             dplyr::transmute(Sample, Study.Group =  case_when(Study.Group %in% case.names ~ "0", 
                                                               .default = "1")) 
    )
  })
  
  # Choose original compounds
  datasets.cmpd <- c()
  covariates_adjust_lst <- list()
  for(d in datasets){
    if(length(data.for.lm[[d]]$cmpd.name) != length(unique(data.for.lm[[d]]$cmpd.name))){
      otu_data <- matrix(0, nrow = nrow(data.for.lm[[d]]$rel.abd), ncol = length(Genera.id),
                         dimnames = list(rownames(data.for.lm[[d]]$rel.abd), Genera.id))
      datasets.cmpd <- c(datasets.cmpd, d)
      cmpd.name <- names(which(table(data.for.lm[[d]]$cmpd.name) > 1))
      cmpd.dat <- data.for.lm[[d]]$rel.cmpd[,data.for.lm[[d]]$cmpd.name %in% cmpd.name]
      otu_data[,colnames(data.for.lm[[d]]$rel.abd)] <- data.for.lm[[d]]$rel.abd %>% as.matrix()
      sample.kep <- intersect(rownames(otu_data), rownames(cmpd.dat))
      cmpd.lst <- c()
      for(l in cmpd.name){
        cmpd.lst <- rbind(cmpd.lst, cbind(colnames(data.for.lm[[d]]$rel.cmpd)[data.for.lm[[d]]$cmpd.name == l],l))
      }
      colnames(cmpd.lst) <- c("Orig.Compound", "Compound")
      cmpd.lst <- data.frame(cmpd.lst)
      covariates_adjust_lst[[d]] <- covariates.disease[[d]][match(sample.kep,covariates.disease[[d]]$Sample),] %>% 
        data.frame(row.names = NULL) %>% tibble::column_to_rownames(var = "Sample")
      otu_data_single <- list(otu_data[sample.kep,])
      names(otu_data_single) <- d
      if(d != "POYET_BIO_ML_2019"){
        null.obj.single <- melody.null.model(rel.abd = otu_data_single, covariate.adjust = covariates_adjust_lst[d])
      }else{
        null.obj.single <- melody.null.model(rel.abd = otu_data_single)
      }
      cluster.tmp <- list(cluster[[d]]$Subject)
      names(cluster.tmp) <- d
      names(cluster.tmp[[1]]) <- rownames(cluster[[d]])
      covariate.interest.tmp <- list(cmpd.dat[sample.kep,])
      names(covariate.interest.tmp) <- d

      summary.stat.study.single <- melody.get.summary(null.obj = null.obj.single,
                                                      covariate.interest = covariate.interest.tmp,
                                                      cluster = cluster.tmp)

      Melody.scan.model <- melody.meta.summary(summary.stats = summary.stat.study.single,
                                               verbose = TRUE)

      cmpd_sty_scan <- data.frame(selected_num = unlist(lapply(Melody.scan.model, function(d){sum(d$coef!=0)}))) %>%
        tibble::rownames_to_column(var = "Orig_Compound") %>% 
        dplyr::mutate(`Orig.Compound` = sort(cmpd.lst$Orig.Compound)) %>%
        dplyr::left_join(cmpd.lst, by = "Orig.Compound")
      
      save(cmpd_sty_scan, file = paste0("./Metabolites/Metabolite_select/Compound_", d, ".Rdata"))
    }
  }
  
  ## Process data
  for(d in datasets){
    if(d %in% datasets.cmpd){
      load(paste0("./Metabolites/Metabolite_select/Compound_", d, ".Rdata"))
      remove.cmpd <- setdiff(cmpd_sty_scan$Orig.Compound, cmpd_sty_scan %>% 
                               dplyr::filter(selected_num == max(selected_num), .by = Compound) %>% 
                               dplyr::filter(row_number() == 1, .by = Compound) %>% 
                               pull(Orig.Compound))
      keep.cmpd <- !(colnames(data.for.lm[[d]]$rel.cmpd) %in% remove.cmpd)
      cmpd.dat <- data.for.lm[[d]]$rel.cmpd[,keep.cmpd]
      data.for.lm[[d]]$cmpd.name <- data.for.lm[[d]]$cmpd.name[keep.cmpd]
      colnames(cmpd.dat) <- data.for.lm[[d]]$cmpd.name 
      data.for.lm[[d]]$rel.cmpd <- cmpd.dat
    }else{
      colnames(data.for.lm[[d]]$rel.cmpd) <- data.for.lm[[d]]$cmpd.name 
    }
  }
  
  otu_data_lst <- list()
  cmpd_data_lst <- list()
  cluster_data_lst <- list()
  covariates_adjust_lst <- list()
  for(d in datasets){
    sample.kep <- intersect(rownames(data.for.lm[[d]]$rel.abd), rownames(data.for.lm[[d]]$rel.cmpd))
    covariates_adjust_lst[[d]] <- covariates.disease[[d]][match(sample.kep, covariates.disease[[d]]$Sample),] %>%
      data.frame(row.names = NULL) %>%
      tibble::column_to_rownames(var = "Sample")
    otu_data_lst[[d]] <- data.for.lm[[d]]$rel.abd %>% dplyr::slice(match(sample.kep, rownames(data.for.lm[[d]]$rel.abd)))
    cmpd_data_lst[[d]] <- data.for.lm[[d]]$rel.cmpd %>% dplyr::slice(match(sample.kep, rownames(data.for.lm[[d]]$rel.cmpd))) 
    cluster_data_lst[[d]] <- (cluster[[d]])[sample.kep,]
    names(cluster_data_lst[[d]]) <- sample.kep
  }
  covariates_adjust_lst$POYET_BIO_ML_2019 <- NULL
  
  # Melody meta-analysis
  null.obj <- melody.null.model(rel.abd = otu_data_lst, covariate.adjust = covariates_adjust_lst)
  
  save(null.obj, file = "./Metabolites/Results/Melody.null.model.Rdata")
  
  system.time(
    ## 788.680 s (1500s)
    summary.stat.study.all <- melody.get.summary(null.obj = null.obj, 
                                                 covariate.interest = cmpd_data_lst,
                                                 cluster = cluster_data_lst, 
                                                 verbose = TRUE)
  )
  save(summary.stat.study.all, file = "./Metabolites/Results/summary.stat.study.all.Rdata")
  
  system.time(
    ## 1258.877s / 5018.976 
    Melody.meta.model <- melody.meta.summary(summary.stats = summary.stat.study.all, 
                                             tune.path = "sequence", 
                                             tune.size.sequence = 0:65,
                                             output.best.one = FALSE,
                                             verbose = FALSE)
  )
  
  save(Melody.meta.model, file = "./Metabolites/Results/Melody.meta.all.Rdata")
  
  

  # MMUPHin meta-analysis
  source("./utility/mmuphin.R")
  MMUPHin.meta.pval <- data.frame(feature = Genera.id)
  MMUPHin.meta.coef <- data.frame(feature = Genera.id)
  for(cmpds.id in Metabolite.ids){
    otu_data_lst <- NULL
    cmpd_data_lst <- NULL
    cluster_data_lst <- NULL
    covariates_adjust_lst <- NULL
    ref.studies <- c()
    for(d in datasets){
      if(cmpds.id %in% colnames(data.for.lm[[d]]$rel.cmpd)){
        sample.kep <- intersect(rownames(data.for.lm[[d]]$rel.abd), 
                                rownames(data.for.lm[[d]]$rel.cmpd)[!is.na(data.for.lm[[d]]$rel.cmpd[,cmpds.id])])
        covariates_adjust_lst <- rbind(covariates_adjust_lst, covariates.disease[[d]][match(sample.kep, covariates.disease[[d]]$Sample),] %>% 
                                         data.frame(row.names = NULL) %>%
                                         tibble::column_to_rownames(var = "Sample"))
        
        tmp.otu <- matrix(0, nrow = length(sample.kep), ncol = length(Genera.id), 
                          dimnames = list(sample.kep, Genera.id))
        inter.taxa <- colnames(null.obj[[d]]$res)
        tmp.otu[sample.kep, inter.taxa] <- as.matrix(data.for.lm[[d]]$rel.abd[sample.kep,inter.taxa])
        otu_data_lst <- rbind(otu_data_lst, tmp.otu)
        tmp.cmpd <- matrix(0, nrow = length(sample.kep), ncol = 1, 
                           dimnames = list(sample.kep, cmpds.id))
        tmp.cmpd[sample.kep, cmpds.id] <- as.matrix(data.for.lm[[d]]$rel.cmpd[sample.kep,cmpds.id])
        cmpd_data_lst <- rbind(cmpd_data_lst, tmp.cmpd)
        cluster_data_lst <- rbind(cluster_data_lst, cluster[[d]])
        ref.studies <- c(ref.studies, rep(d, nrow(tmp.otu)))
      }
    }
    cluster_data_lst <- cluster_data_lst[rownames(otu_data_lst),]
    meta.data <- cbind(cluster_data_lst, covariates_adjust_lst, ref.studies, cmpd_data_lst) 
    colnames(meta.data) <- c("Subject", "Study.Group", "Study" ,cmpds.id)
    meta.data$Study <- factor(meta.data$Study)
    
    if(length(unique(meta.data$Subject)) == length(meta.data$Subject)){
      MMUPHin.model <- fit_metaAnalysis(feature.abd = t(round(otu_data_lst)), #[,!(gsub(".*;g__","g__", colnames(otu_data_lst)) %in% rm.genus)])),
                                        data = meta.data,
                                        test_variable = cmpds.id,
                                        batch_variable = "Study",
                                        covariates = "Study.Group",
                                        covariates.random = NULL,
                                        moderator_variables = NULL)
    }else{
      MMUPHin.model <- fit_metaAnalysis(feature.abd = t(round(otu_data_lst)), #[,!(gsub(".*;g__","g__", colnames(otu_data_lst)) %in% rm.genus)])),
                                        data = meta.data,
                                        test_variable = cmpds.id,
                                        batch_variable = "Study",
                                        covariates = "Study.Group",
                                        covariates.random = "Subject",
                                        moderator_variables = NULL)
    }
    
    sign.est <- matrix(NA, nrow = length(Genera.id), ncol = length(names(MMUPHin.model$maaslin_fits)), 
                       dimnames = list(Genera.id, names(MMUPHin.model$maaslin_fits)))
    for(d in names(MMUPHin.model$maaslin_fits)){
      sign.est[MMUPHin.model$maaslin_fits[[d]]$feature,d] <- 
        sign(MMUPHin.model$maaslin_fits[[d]]$coef) == sign(MMUPHin.model$meta_fits[MMUPHin.model$maaslin_fits[[d]]$feature,]$coef)
    }
    sign.ratio <- rowMeans(sign.est, na.rm = TRUE)
    
    MMUPHin.meta.pval <- MMUPHin.meta.pval %>% 
      dplyr::left_join(MMUPHin.model$meta_fits %>% 
                         dplyr::transmute(feature, pval) %>% 
                         dplyr::rename_with(~cmpds.id, pval), by = "feature")
    
    MMUPHin.meta.coef <- MMUPHin.meta.coef %>% 
      dplyr::left_join(MMUPHin.model$meta_fits %>% 
                         dplyr::transmute(feature, coef) %>%
                         dplyr::rename_with(~cmpds.id, coef), by = "feature")
    
  }
  
  save(MMUPHin.meta.coef, MMUPHin.meta.pval, file = "./Metabolites/Results/MMUPHin.meta.all.Rdata")
  
  
  