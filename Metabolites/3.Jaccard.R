

  rm(list = ls())
  ## Package
  library("miMeta")
  library("tidyverse")
  library("ggplot2")
  
  ## Read data
  load("./Metabolites/Processed_data/processed_data.Rdata")
  
  ## process data
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
  
  #################################### Choose original compounds ####################################
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
        cmpd.lst <- rbind(cmpd.lst,
                          cbind(colnames(data.for.lm[[d]]$rel.cmpd)[data.for.lm[[d]]$cmpd.name == l],
                                l))
      }
      colnames(cmpd.lst) <- c("Orig.Compound", "Compound")
      cmpd.lst <- data.frame(cmpd.lst)
      covariates_adjust_lst[[d]] <- covariates.disease[[d]][match(sample.kep,covariates.disease[[d]]$Sample),] %>%
        data.frame(row.names = NULL) %>% tibble::column_to_rownames(var = "Sample")
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
  
  ## Reform data
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
  
  ## Compute Jaccard Index before an after remove top genera
  load("./Metabolites/Results/Melody.meta.all.Rdata")
  load("./Metabolites/Results/MMUPHin.meta.all.Rdata")
  feature.ids <- rownames(Melody.meta.model$HMDB0000008$coef)
  Melody_coef_mat <- matrix(0, nrow = length(feature.ids), ncol = length(Melody.meta.model),
                            dimnames = list(feature.ids, names(Melody.meta.model)))
  for(ids in names(Melody.meta.model)){
    Melody_coef_mat[rownames(Melody.meta.model[[ids]]$coef), ids] <- Melody.meta.model[[ids]]$coef[,which.min(Melody.meta.model[[ids]]$dev + Melody.meta.model[[ids]]$ic)]
  }
  Melody_coef_mat[is.na(Melody_coef_mat)] <- 0
  MMUPHin_coef_mat <- MMUPHin.meta.coef %>% tibble::column_to_rownames(var = "feature")
  MMUPHin_pval_mat <- MMUPHin.meta.pval %>% tibble::column_to_rownames(var = "feature")
  MMUPHin_qval_mat <- apply(MMUPHin_pval_mat, 2, p.adjust, method = "fdr")
  MMUPHin_coef_mat[MMUPHin_qval_mat > 0.05] <- 0
  MMUPHin_coef_mat[is.na(MMUPHin_coef_mat)] <- 0
  
  ## Pick Metabolites
  rm.cp.lst <- c()
  selected.len <- c()
  for(l in colnames(MMUPHin_coef_mat)){
    if(sum(MMUPHin_coef_mat[,l]>0) >= 5 & sum(Melody_coef_mat[,l]>0) >= 5 & 
       sum(MMUPHin_coef_mat[,l]!=0) > 5 & sum(Melody_coef_mat[,l]!=0) > 5){
      rm.cp.lst <- c(rm.cp.lst, l)
    }
  }
  
  Melody.bf <- data.frame(feature = Genera.id)
  Melody.af.match <- data.frame(feature = Genera.id)
  Melody.af.best <- data.frame(feature = Genera.id)
  for(l in rm.cp.lst){
    ## before genus removal
    otu_data_lst_sub <- list()
    cmpd_data_lst_sub <- list()
    for(d in names(otu_data_lst)){
      if(l %in% colnames(cmpd_data_lst[[d]])){
        cmpd_data_lst_sub[[d]] <- cmpd_data_lst[[d]] %>% select(all_of(l))
        otu_data_lst_sub[[d]] <- otu_data_lst[[d]]
      }
    }
    null.obj <- melody.null.model(rel.abd = otu_data_lst_sub, covariate.adjust = covariates_adjust_lst)
    
    summary.stat.study <- melody.get.summary(null.obj = null.obj, 
                                             covariate.interest = cmpd_data_lst_sub,
                                             cluster = cluster_data_lst)
    
    Melody.model <- melody.meta.summary(summary.stats = summary.stat.study, 
                                        tune.path = "sequence",
                                        tune.size.sequence = 0:65,
                                        output.best.one = FALSE,
                                        verbose = TRUE)
    
    
    Melody.bf <- Melody.bf %>% 
      left_join(data.frame(coef = Melody.model[[l]]$coef[,which.min(Melody.model[[l]]$dev + Melody.model[[l]]$ic)]) %>% 
                  tibble::rownames_to_column("feature") %>%
                  dplyr::rename_with(~l, coef), by = "feature")
    
    top5.genera <- Melody.bf$feature[rank(-Melody.bf[[l]]) <= 5]
    
    rm(list = c("null.obj", "summary.stat.study", "Melody.model"))
    
    ## after genus removal
    otu_data_lst_sub <- list()
    cmpd_data_lst_sub <- list()
    for(d in names(otu_data_lst)){
      if(l %in% colnames(cmpd_data_lst[[d]])){
        cmpd_data_lst_sub[[d]] <- cmpd_data_lst[[d]] %>% select(all_of(l))
        otu_data_lst_sub[[d]] <- (otu_data_lst[[d]])[,!(colnames(otu_data_lst[[d]]) %in% top5.genera)]
      }
    }
    null.obj <- melody.null.model(rel.abd = otu_data_lst_sub, covariate.adjust = covariates_adjust_lst)
    
    summary.stat.study <- melody.get.summary(null.obj = null.obj, 
                                             covariate.interest = cmpd_data_lst_sub,
                                             cluster = cluster_data_lst)
    
    Melody.model <- melody.meta.summary(summary.stats = summary.stat.study, 
                                        tune.path = "sequence",
                                        tune.size.sequence = 0:60,
                                        output.best.one = FALSE,
                                        verbose = TRUE)
    
    Melody.af.best <- Melody.af.best %>% 
      left_join(data.frame(coef = Melody.model[[l]]$coef[,which.min(Melody.model[[l]]$dev + Melody.model[[l]]$ic)]) %>% 
                  tibble::rownames_to_column("feature") %>%
                  dplyr::rename_with(~l, coef), by = "feature")
    
    Melody.af.match <- Melody.af.match %>% 
      left_join(data.frame(coef = Melody.model[[l]]$coef[,sum(Melody.bf[[l]]!=0) - length(top5.genera) + 1]) %>% 
                  tibble::rownames_to_column("feature") %>%
                  dplyr::rename_with(~l, coef), by = "feature")
    
    rm(list = c("null.obj", "summary.stat.study", "Melody.model"))
  }
  
  save(Melody.bf, Melody.af.best, Melody.af.match, file = "./Metabolites/Results/Melody.Jaccard.Rdata")
  
  ## MMUPHin default model
  source("./utility/mmuphin.R")
  load("./Metabolites/Results/Melody.null.model.Rdata")
  
  MMUPHin.bf.pval <- data.frame(feature = Genera.id)
  MMUPHin.bf.coef <- data.frame(feature = Genera.id)
  MMUPHin.af.pval <- data.frame(feature = Genera.id)
  MMUPHin.af.coef <- data.frame(feature = Genera.id)
  for(l in rm.cp.lst){
    otu_data_lst <- NULL
    cmpd_data_lst <- NULL
    cluster_data_lst <- NULL
    covariates_adjust_lst <- NULL
    ref.studies <- c()
    for(d in datasets){
      if(l %in% colnames(data.for.lm[[d]]$rel.cmpd)){
        sample.kep <- intersect(rownames(data.for.lm[[d]]$rel.abd), 
                                rownames(data.for.lm[[d]]$rel.cmpd)[!is.na(data.for.lm[[d]]$rel.cmpd[,l])])
        covariates_adjust_lst <- rbind(covariates_adjust_lst, covariates.disease[[d]][match(sample.kep, covariates.disease[[d]]$Sample),] %>% 
                                         data.frame(row.names = NULL) %>%
                                         tibble::column_to_rownames(var = "Sample"))
        
        tmp.otu <- matrix(0, nrow = length(sample.kep), ncol = length(Genera.id), 
                          dimnames = list(sample.kep, Genera.id))
        inter.taxa <- colnames(null.obj[[d]]$res)
        tmp.otu[sample.kep, inter.taxa] <- as.matrix(data.for.lm[[d]]$rel.abd[sample.kep,inter.taxa])
        otu_data_lst <- rbind(otu_data_lst, tmp.otu)
        tmp.cmpd <- matrix(0, nrow = length(sample.kep), ncol = 1, 
                           dimnames = list(sample.kep, l))
        tmp.cmpd[sample.kep, l] <- as.matrix(data.for.lm[[d]]$rel.cmpd[sample.kep,l])
        cmpd_data_lst <- rbind(cmpd_data_lst, tmp.cmpd)
        cluster_data_lst <- rbind(cluster_data_lst, cluster[[d]])
        ref.studies <- c(ref.studies, rep(d, nrow(tmp.otu)))
      }
    }
    cluster_data_lst <- cluster_data_lst[rownames(otu_data_lst),]
    meta.data <- cbind(cluster_data_lst, covariates_adjust_lst, ref.studies, cmpd_data_lst) 
    colnames(meta.data) <- c("Subject", "Study.Group", "Study" ,l)
    meta.data$Study <- factor(meta.data$Study)
    
    if(length(unique(meta.data$Subject)) == length(meta.data$Subject)){
      MMUPHin.model <- fit_metaAnalysis(feature.abd = t(round(otu_data_lst)),
                                        data = meta.data,
                                        test_variable = l,
                                        batch_variable = "Study",
                                        covariates = "Study.Group",
                                        covariates.random = NULL,
                                        moderator_variables = NULL)
    }else{
      MMUPHin.model <- fit_metaAnalysis(feature.abd = t(round(otu_data_lst)),
                                        data = meta.data,
                                        test_variable = l,
                                        batch_variable = "Study",
                                        covariates = "Study.Group",
                                        covariates.random = "Subject",
                                        moderator_variables = NULL)
    }
    
    
    MMUPHin.bf.pval <- MMUPHin.bf.pval %>% dplyr::left_join(MMUPHin.model$meta_fits %>% 
                                                              dplyr::transmute(feature, pval) %>% 
                                                              dplyr::rename_with(~l, pval), by = "feature")
    
    MMUPHin.bf.coef <- MMUPHin.bf.coef %>% dplyr::left_join(MMUPHin.model$meta_fits %>% 
                                                              dplyr::transmute(feature, coef) %>%
                                                              dplyr::rename_with(~l, coef), by = "feature")
    
    MMUPHin.tmp.coef <- MMUPHin.bf.coef[[l]]
    MMUPHin.tmp.coef[p.adjust(MMUPHin.bf.pval[[l]], method = "BH") > 0.05] <- 0
    top5.genera <- MMUPHin.bf.coef$feature[rank(-MMUPHin.tmp.coef) <= 5]
    
    if(length(unique(meta.data$Subject)) == length(meta.data$Subject)){
      MMUPHin.model <- fit_metaAnalysis(feature.abd = t(round(otu_data_lst[,!(colnames(otu_data_lst) %in% top5.genera)])),
                                        data = meta.data,
                                        test_variable = l,
                                        batch_variable = "Study",
                                        covariates = "Study.Group",
                                        covariates.random = NULL,
                                        moderator_variables = NULL)
    }else{
      MMUPHin.model <- fit_metaAnalysis(feature.abd = t(round(otu_data_lst[,!(colnames(otu_data_lst) %in% top5.genera)])),
                                        data = meta.data,
                                        test_variable = l,
                                        batch_variable = "Study",
                                        covariates = "Study.Group",
                                        covariates.random = "Subject",
                                        moderator_variables = NULL)
    }
    
    
    MMUPHin.af.pval <- MMUPHin.af.pval %>% dplyr::left_join(MMUPHin.model$meta_fits %>% 
                                                              dplyr::transmute(feature, pval) %>% 
                                                              dplyr::rename_with(~l, pval), by = "feature")
    
    MMUPHin.af.coef <- MMUPHin.af.coef %>% dplyr::left_join(MMUPHin.model$meta_fits %>% 
                                                              dplyr::transmute(feature, coef) %>%
                                                              dplyr::rename_with(~l, coef), by = "feature")
    
  }

  save(MMUPHin.af.pval, MMUPHin.af.coef, MMUPHin.bf.pval, MMUPHin.bf.coef, file = "./Metabolites/Results/MMUPHin.Jaccard.Rdata")
  
  