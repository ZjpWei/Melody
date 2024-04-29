# =============================================== #
#    (1) Prepare microbiome-metabolite data       #
# =============================================== #

  library(dplyr)
  library(ggplot2)
  library(logger)
  
  # Notebook settings ----
  future::plan("multisession", workers = 4)
  options(scipen = 999)
  
  ## Load utility scripts
  source("./utility/utils.R")
  
  ## Load all data available in the curated gut microbiome-metabolome data resource:  
  all.data <- load.all.datasets(parent.folder = "./Data/Metabolite_data/processed_data")
  for(i in 1:length(all.data)) assign(names(all.data)[i], all.data[[i]])
  rm(all.data)
  
  ### Analysis settings
  # Limit the analysis to non-infants cohorts only.We allow multiple samples per 
  # subject (in longitudinal studies).  
  
  datasets <- c("ERAWIJANTARI_GASTRIC_CANCER_2020",
                "YACHIDA_CRC_2019",
                "KIM_ADENOMAS_2020",
                "FRANZOSA_IBD_2019",
                "MARS_IBS_2020",
                "iHMP_IBDMDB_2019",
                "WANG_ESRD_2020",
                "POYET_BIO_ML_2019")
  
  # Remove subjects from YACHIDA_CRC_2019 study that are also in ERAWIJANTARI_GASTRIC_CANCER_2020 study:
  metadata$YACHIDA_CRC_2019 <- metadata$YACHIDA_CRC_2019 %>%
    dplyr::filter(! Shared.w.ERAWIJANTARI_2020)
  updated.yachida.sample.list <- metadata$YACHIDA_CRC_2019$Sample
  mtb$YACHIDA_CRC_2019 <- mtb$YACHIDA_CRC_2019 %>% filter(Sample %in% updated.yachida.sample.list)
  genera$YACHIDA_CRC_2019 <- genera$YACHIDA_CRC_2019 %>% filter(Sample %in% updated.yachida.sample.list)
  metadata$POYET_BIO_ML_2019$Study.Group <- "NA"
  
  # Genera plot and compounds plots
  genera.dataset.stats <- get.genera.dataset.stats(genera, datasets)
  
  # We further summarize these stats to get basic statistics 
  #  at the genus level (average over datasets)
  genera.stats <- genera.dataset.stats %>%
    group_by(Taxon) %>%
    summarise(Taxon.Overall.Mean.Abundance = 
                weighted.mean(x = Taxon.Mean.Abundance, w = Dataset.N),
              Taxon.Overall.Perc.Non.Zero = 
                weighted.mean(x = Taxon.Perc.of.Non.Zeros, w = Dataset.N),
              N.Datasets.Including.Taxon = n(),
              .groups = "drop") %>%
    mutate(Genus.Only = gsub(".*\\;g__","g__", Taxon))
  
  # Print/plot statistics
  message(paste(nrow(genera.stats), "unique genera (or higher-level clades) were found across all datasets"))
  
  tmp <- genera.stats %>%
    group_by(N.Datasets.Including.Taxon) %>%
    summarise(N = n(), .groups = "drop") %>%
    arrange(-N.Datasets.Including.Taxon) %>%
    mutate(cum.N = cumsum(N))
  
  label1.x <- unlist(tmp[tmp$N.Datasets.Including.Taxon == 1, "cum.N"])
  label2.x <- unlist(tmp[tmp$N.Datasets.Including.Taxon == 10, "cum.N"])
  label1 <- paste(label1.x,"genera in total")
  label2 <- paste(label2.x,"genera are present\nin at least 10 datasets")
  label.y <- max(tmp$cum.N)
  
  # pdf("./figure_genera_bar.pdf", width = 5.25, height = 4.8, bg = "white")
  
  ggplot(tmp, aes(x = N.Datasets.Including.Taxon, y = cum.N)) +
    geom_bar(stat="identity", color='black', 
             fill='#8CBA80') +
    theme_classic() +
    ylab("No. of genera") +
    xlab("No. of datasets") +
    scale_x_continuous(breaks=2*1:4, limits = c(NA,8.5)) +
    annotate(geom = "curve", 
             x = c(3,10),
             y = c(label1.x,label2.x)+800, 
             xend = c(3.5,10.6), 
             yend = c(0.99*label.y, 0.32*label.y), 
             color = "#7E7F9A",
             curvature = -.3, 
             arrow = arrow(length = unit(2, "mm"))) +
    annotate("text", 
             x = c(3.7, 10.8), 
             y = c(0.97*label.y, 0.3*label.y)+100, 
             label = c(label1, label2), 
             color = "#272838", 
             size = 3.45, hjust = 0) +
    annotate(geom = "point",
             x = c(3,10), 
             y = c(label1.x,label2.x)+800,
             color = "#7E7F9A", size = 2) +
    geom_text(aes(label=cum.N), vjust = -0.5, size = 2.6) +
    scale_y_continuous(limits = c(0, max(tmp$cum.N) + 7)) +
    theme(axis.title = element_text(size = 11)) 
  
  # dev.off()
  
  rm(label1.x, label2.x, label1, label2, label.y, tmp)
  
  metabolites.per.dataset <- 
    get.metab.dataset.stats(mtb.map, datasets)
  
  # Count for each metabolite how many datasets have it
  metabolites.stats <- metabolites.per.dataset %>%
    group_by(Type, Compound) %>%
    summarise(N = n(), N.Datasets.Including.Compound = n_distinct(Dataset), .groups = "drop")
  
  # Print statistics 
  print(paste(n_distinct(metabolites.stats$Compound[metabolites.stats$Type == "HMDB"]),
              "unique HMDB compound IDs were found across all datasets"))
  print(paste(n_distinct(metabolites.stats$Compound[metabolites.stats$Type == "KEGG"]),
              "unique KEGG compound IDs were found across all datasets"))
  
  tmp <- metabolites.stats %>%
    group_by(Type, N.Datasets.Including.Compound) %>%
    summarise(N = n(), .groups = "drop") %>%
    group_by(Type) %>%
    arrange(-N.Datasets.Including.Compound) %>%
    mutate(cum.N = cumsum(N))
  
  # Focus on HMDB
  tmp <- tmp %>% filter(Type == "HMDB")
  
  label1.x <- unlist(tmp[tmp$N.Datasets.Including.Compound == 1, "cum.N"])
  label2.x <- unlist(tmp[tmp$N.Datasets.Including.Compound == 10, "cum.N"])
  label1 <- paste(label1.x,"compounds \n (HMDB-annotated) in total")
  label2 <- paste(label2.x,"compounds (HMDB)\nare present in at least\n10 datasets")
  label.y <- max(tmp$cum.N)
  
  # pdf("./figure_cmpd_bar.pdf", width = 5.25, height = 4.8, bg = "white")
  
  ggplot(tmp, aes(x = N.Datasets.Including.Compound, y = cum.N)) +
    geom_bar(stat="identity", color='black', 
             fill='#FDAE6B') +
    theme_classic() +
    ylab("No. of compounds") +
    xlab("No. of datasets") +
    scale_x_continuous(breaks=2*1:4, limits = c(NA,8.5)) +
    annotate(geom = "curve", 
             x = c(3,10), 
             y = c(label1.x,label2.x)+100, 
             xend = c(3,10)+0.8, 
             yend = c(0.7*label.y, 0.5*label.y), 
             color = "#7E7F9A",
             curvature = -.3, 
             arrow = arrow(length = unit(2, "mm"))) +
    annotate("text", 
             x = c(4, 11), 
             y = c(0.7*label.y, 0.5*label.y)-30, 
             label = c(label1, label2), 
             color = "#272838", 
             size = 3.45, hjust = 0) +
    annotate(geom = "point",
             x = c(3,10), y = c(label1.x,label2.x)+100,
             color = "#7E7F9A", size = 2) +
    geom_text(aes(label=cum.N), vjust = -0.5, size = 2.6) +
    scale_y_continuous(limits = c(0, max(tmp$cum.N) + 15)) +
    theme(axis.title = element_text(size = 11))
  
  # dev.off()
  
  rm(label1.x, label2.x, label1, label2, label.y, tmp)
  
  # Remove compound which has less than 0.1 prevalence
  for(d in datasets){
    tmp.mbp <- mtb[[d]] %>% tibble::column_to_rownames("Sample")
    ## convert NA to 0 
    tmp.mbp <- tmp.mbp %>% dplyr::mutate_all(~replace(., is.na(.), 0))
    ## 0.1 prevalence filter
    filter.cmpd <- names(which(colMeans(tmp.mbp !=0, na.rm = TRUE) >= 0.1))
    tmp.mbp <- tmp.mbp %>% tibble::rownames_to_column("Sample")
    mtb[[d]] <- tmp.mbp[,c("Sample", filter.cmpd)]
    mtb.map[[d]] <- (mtb.map[[d]])[mtb.map[[d]]$Compound %in% filter.cmpd,]
  }
  
  ## Get genus-metabolite pairs 
  # Here we prepare a list of genus-metabolite pairs that appear in at least 2 datasets. 
  # We also require that the genera are not rare (see definition below) and that metabolites 
  # have an HMDB identification and are not constant over samples.  
  # We start by marking for each genus and each metabolite which datasets they are in:
  # Metabolite-dataset statistics
  metabolites.per.dataset <- get.metab.dataset.stats(mtb.map, datasets)
  
  # Genera-dataset statistics
  genera.dataset.stats <- 
    get.genera.dataset.stats(genera, datasets) %>%
    # Add averaged statistics (over datasets)
    group_by(Taxon) %>%
    mutate(Averaged.Taxon.Mean.Abundance = 
             weighted.mean(Taxon.Mean.Abundance, Dataset.N),
           Averaged.Taxon.Perc.of.Non.Zeros = 
             weighted.mean(Taxon.Perc.of.Non.Zeros, Dataset.N),
           N.Datasets = n_distinct(Dataset))
  
  # Discard rare genera (defined here as <25% non-zero values or average abundance 
  # <0.1% over all datasets in this analysis):  
  genera.dataset.stats <- genera.dataset.stats %>%
    filter(Averaged.Taxon.Perc.of.Non.Zeros >= 25) %>%
    filter(Averaged.Taxon.Mean.Abundance >= 0.001)
  
  # We additionally discard genera from individual datasets if they 
  #  are mostly zero's there. See for example:
  #  View(genera.dataset.stats %>% filter(grepl("g__Clostridioides",Taxon)))
  genera.dataset.stats <- genera.dataset.stats %>%
    filter(Taxon.Perc.of.Non.Zeros >= 10) 
  
  # And lastly discard ambiguous/unidentified genera
  genera.dataset.stats <- genera.dataset.stats %>%
    filter(! grepl("g__$", Taxon))
  
  # Discard metabolites with no HMDB annotation or metabolites with constant values:  
  metabolites.per.dataset <- metabolites.per.dataset %>%
    filter(Type == "HMDB") %>%
    select(-Type)
  
  # Also remove metabolites with constant values across cohort
  is.constant <- apply(metabolites.per.dataset, MARGIN = 1, function(r) {
    # Get vector of values of a metabolite
    tmp <- mtb[[unname(r["Dataset"])]][,unname(r["Orig.Compound"])]
    
    # Return true if constant
    return(var(tmp, na.rm = TRUE) == 0)
  })
  metabolites.per.dataset <- metabolites.per.dataset[!is.constant,]
  
  # Retrieve a list of genus-metabolite *pairs* that appear in at least 2 datasets. Save the list in the `common.pairs` table.  
  # Note: some metabolites may appear more than once in a dataset (for example in the case of low-confidence 
  # in annotation or multiple MS runs). We deal with these cases later on
  
  common.pairs <- inner_join(genera.dataset.stats, metabolites.per.dataset, by = "Dataset") %>% 
    relocate(Dataset, Dataset.N) %>% group_by(Taxon, Compound) %>% filter(n_distinct(Dataset) >= 2) %>% 
    mutate(Pair = paste(Compound, gsub(".*;f__","f__",Taxon), sep = "~"))
  
  # Print statistics
  paste(n_distinct(common.pairs$Pair), 
        "unique genus-metabolite pairs will be analyzed")
  paste("These include", n_distinct(common.pairs$Compound), 
        "metabolites and", n_distinct(common.pairs$Taxon), "genera")
  
  metadata.fields <- c("Sample", "Age", "Gender", "Subject", "Study.Group", "BMI")
  
  # Ranked-based inverse normal transformation for metabolites
  mtb.rank <- lapply(mtb, function(d){
    a <- apply(d %>% tibble::column_to_rownames("Sample"), 2,
               function(X){
                 d.rank <- c()
                 for(l in 1:length(X)){
                   d.rank <- c(d.rank, (1 + sum(X[l] > X)))
                 }
                 d.inverse.rank <- qnorm((d.rank - 0.5)/sum(!is.na(d.rank)))
                 return(d.inverse.rank)
               })
    rownames(a) <- d$Sample
    a <- data.frame(a) %>% tibble::rownames_to_column("Sample")
    colnames(a) <- colnames(d)
    return(a)
  })
  
  data.for.lm <- lapply(datasets, function(d) {
    log_info(sprintf("Preparing data for %s", d))
    
    # Get relevant genera
    relevant.genera <- common.pairs %>% filter(Dataset == d) %>% pull(Taxon) %>% unique()
    relevant.genera <- genera.counts[[d]] %>%
      select("Sample", any_of(relevant.genera)) %>% tibble::column_to_rownames("Sample")
    
    # Get relevant metabolites
    relevant.cmpd <- common.pairs %>% filter(Dataset == d) %>% pull(Orig.Compound) %>% unique()
    relevant.cmpd <- mtb.rank[[d]] %>% select("Sample", any_of(relevant.cmpd))
    
    relevant.cmpd <- relevant.cmpd %>% tibble::column_to_rownames("Sample")
    
    # Combine all variables in one table
    tmp <- list(rel.abd = relevant.genera, rel.cmpd = relevant.cmpd, 
                cmpd.name = data.frame(Compound = colnames(relevant.cmpd)) %>% 
                  dplyr::left_join(mtb.map[[d]], by = "Compound") %>%
                  dplyr::pull(HMDB))
    return(tmp)
  })
  names(data.for.lm) <- datasets
  
  save(metadata, data.for.lm, datasets, common.pairs, file = "./Metabolites/Processed_data/processed_data.Rdata")
  
  
