require(readr)

# ----------------------------------------------------------------
# Utility functions for data processing and data analysis
# ----------------------------------------------------------------

# Load all datasets within @parent.folder.
# @parent.folder is expected to contain a child folder per dataset,
#  and the folders name is taken as dataset name.
# The returned object is a list of data types (mtb, mtb.map, genera, 
#  species, metadata) where each data type is a list of tables 
#  named by dataset name
load.all.datasets <- function(parent.folder = "processed_data") {
  # Get all processed datasets
  data.dirs <- list.dirs(file.path(parent.folder))[-1]
  
  # Initialize table lists
  all.data <- list()
  all.data$data.dirs <- data.dirs
  all.data$metadata <- list()
  all.data$mtb <- list()
  all.data$mtb.map <- list()
  all.data$genera <- list()
  all.data$species <- list()
  all.data$genera.counts <- list()
  all.data$species.counts <- list()
  
  for (x in data.dirs) {
    # Create a temporary environment to hold all processed tables
    tmp.env <- new.env()
    dataset.name <- basename(x)
    
    # Load and save tables
    load(file.path(x, ".RData"), tmp.env)
    all.data$mtb[[dataset.name]] <- get('mtb', tmp.env) 
    all.data$mtb.map[[dataset.name]] <- get('mtb.map', tmp.env) 
    all.data$genera[[dataset.name]] <- get('genera', tmp.env) 
    all.data$metadata[[dataset.name]] <- get('metadata', tmp.env)
    if ("species" %in% ls(tmp.env)) all.data$species[[dataset.name]] <- get('species', tmp.env) 
    if ("genera.counts" %in% ls(tmp.env)) all.data$genera.counts[[dataset.name]] <- get('genera.counts', tmp.env) 
    if ("species.counts" %in% ls(tmp.env)) all.data$species.counts[[dataset.name]] <- get('species.counts', tmp.env) 
    
    # Clean up
    rm(tmp.env)
  }
  
  message("Datasets loaded successfully")
  return(all.data)
}

get.genera.dataset.stats <- function(genera, datasets) {
  # Initialize table for genera statistics 
  genera.dataset.stats <- 
    data.frame(Taxon = character(0),
               Dataset = character(0),
               Dataset.N = integer(0),
               Taxon.Mean.Abundance = numeric(0),
               Taxon.Var.Abundance = numeric(0),
               Taxon.Perc.of.Non.Zeros = numeric(0))
  
  # We iterate over each dataset and record a few basic stats 
  #  per each genus in each dataset
  for (dataset in datasets) {
    tmp <- genera[[dataset]] %>% select(-Sample)
    tmp.genera <- colnames(tmp)
    tmp.means <- unname(apply(tmp, MARGIN = 2, mean))
    tmp.vars <- unname(apply(tmp, MARGIN = 2, var))
    tmp.non.zero.perc <- unname(apply(tmp, MARGIN = 2, 
                                      function(v) {100*sum(v>0)/length(v)}))
    genera.dataset.stats <- 
      bind_rows(genera.dataset.stats, 
                data.frame(Taxon = tmp.genera, 
                           Dataset = dataset,
                           Dataset.N = nrow(tmp),
                           Taxon.Mean.Abundance = tmp.means,
                           Taxon.Var.Abundance = tmp.vars,
                           Taxon.Perc.of.Non.Zeros = tmp.non.zero.perc))
  }
  
  return(genera.dataset.stats)
}

get.metab.dataset.stats <- function(mtb.map, datasets) {
  metabolites.per.dataset <- data.frame(stringsAsFactors = F)
  
  for (dataset in datasets) {
    tmp <- mtb.map[[dataset]] %>%
      select(Compound, KEGG, HMDB) %>%
      rename(Orig.Compound = Compound) %>%
      tidyr::pivot_longer(cols = c("KEGG","HMDB"), 
                          names_to = "Type", 
                          values_to = "Compound", 
                          values_drop_na = TRUE) %>%
      mutate(Dataset = dataset) 
    metabolites.per.dataset <- bind_rows(metabolites.per.dataset, tmp)
  }
  
  return(metabolites.per.dataset)
}

# Get significant marks given p values, for plotting
get.signif.marks <- function(p.vals) {
  signif.marks <- sapply(p.vals,
                         function(x) {
                           ifelse(x <= 0.001, "***", 
                                  ifelse(x <= 0.01, "**", 
                                         ifelse(x <= 0.05, "*", 
                                                ifelse(x <= 0.1, "^", ""))))
                         })
  signif.marks[is.na(signif.marks)] <- ""
  return(signif.marks)
}

# Add a new predictor variable to an existing formula object
add.var.to.formula <- function(form, new.var) {
  require(formula.tools)
  new.form.str <- paste(as.character(form), "+", new.var)
  return(as.formula(new.form.str))
}
