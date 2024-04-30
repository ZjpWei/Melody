# =============================================== #
#             Supplementary figure 8              #
# =============================================== #

  # Packages ----
  library(ggplot2)
  library(dplyr)
  library(gt)
  library(logger)
  
  # Notebook settings ----
  future::plan("multisession", workers = 4)
  options(scipen = 999)
  
  ## Load utility scripts
  source("./utility/utils.R")
  
  ## Load all data available in the curated gut microbiome-metabolome data resource:  
  all.data <- load.all.datasets("./Data/Metabolite_data/processed_data")
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
    relocate(Dataset, Dataset.N) %>% group_by(Taxon, Compound) %>% filter(n_distinct(Dataset) >= 1) %>% 
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
    
    # For compactness we fetch only relevant genera and metabolites,
    #  included in our pairs of interest.
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
  
  Compounds <- NULL
  Genera <- NULL
  for(d in names(data.for.lm)){
    Compounds <- c(Compounds, unique(data.for.lm[[d]]$cmpd.name))
    Genera <- c(Genera, colnames(data.for.lm[[d]]$rel.abd))
  }
  
  genera.stats <- data.frame(Count = table(Genera)) %>% 
    dplyr::transmute(Compound = Count.Genera, 
                     N.Datasets.Including.Taxon = Count.Freq)
  
  tmp <- genera.stats %>%
    group_by(N.Datasets.Including.Taxon) %>%
    summarise(N = n(), .groups = "drop") %>%
    arrange(-N.Datasets.Including.Taxon) %>%
    mutate(cum.N = cumsum(N)) %>% 
    add_row(N.Datasets.Including.Taxon = 5:1, N = 0, cum.N = 101)
  
  p1 <- ggplot(tmp, aes(x = N.Datasets.Including.Taxon, y = cum.N)) +
    geom_bar(stat="identity", color='black', 
             fill='gray') +
    theme_classic() +
    ylab("Number of genera") +
    xlab("Number of studies") + 
    scale_x_continuous(breaks=1:8, limits = c(NA,8.5)) +
    geom_text(aes(label=cum.N), vjust = -0.5, size = 2.6) +
    scale_y_continuous(limits = c(0, max(tmp$cum.N) + 7)) +
    theme(axis.title = element_text(size = 11)) 
  
  # Clean up 
  rm(tmp)
  
  metabolites.stats <- data.frame(Count = table(Compounds)) %>% 
    dplyr::transmute(Compound = Count.Compounds, 
                     N.Datasets.Including.Compound = Count.Freq)
  
  tmp <- metabolites.stats %>%
    group_by(N.Datasets.Including.Compound) %>%
    summarise(N = n(), .groups = "drop") %>%
    arrange(-N.Datasets.Including.Compound) %>%
    mutate(cum.N = cumsum(N))
  
  # Focus on HMDB
  
  p2 <- ggplot(tmp, aes(x = N.Datasets.Including.Compound, y = cum.N)) +
    geom_bar(stat="identity", color='black', 
             fill='#FDAE6B') +
    theme_classic() +
    ylab("Number of metabolites") +
    xlab("Number of studies") +
    scale_x_continuous(breaks=1:8, limits = c(NA,8.5)) +
    geom_text(aes(label=cum.N), vjust = -0.5, size = 2.6) +
    scale_y_continuous(limits = c(0, max(tmp$cum.N) + 15)) +
    theme(axis.title = element_text(size = 11))
  
  rm(tmp)
  
  pdf("./figures/FigS8.pdf", width = 6.85, height = 3.33, bg = "white")
  
  ggpubr::ggarrange(p2, p1, nrow = 1, ncol = 2)
  
  dev.off()


