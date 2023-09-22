library("ANCOMBC")

ancombc.fun <- function(feature.table, meta, formula, types = "ancombc") {
  assays <- S4Vectors::SimpleList(counts = as.matrix(feature.table))
  smd <- S4Vectors::DataFrame(meta)
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
  if(types == "ancombc2"){
    model <- ancombc2(data = tse, assay_name = "counts", tax_level = NULL, 
                      fix_formula = formula, rand_formula = NULL,
                      p_adj_method = "fdr", prv_cut = 0, lib_cut = 0,
                      s0_perc = 0, group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                      alpha = 0.01, n_cl = 2, verbose = TRUE, global = FALSE, pairwise = FALSE,
                      dunnet = FALSE, trend = FALSE)
  }else{
    model <- ancombc(data = tse, assay_name = "counts", tax_level = NULL, phyloseq = NULL,
                     formula = formula, p_adj_method = "fdr", prv_cut = 0, lib_cut = 0,
                     group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                     alpha = 0.01, n_cl = 2, verbose = TRUE, global = FALSE)
  }
  return(model)
}
