library("ANCOMBC")

ancombc.fun <- function(feature.table, 
                        meta, 
                        formula, 
                        adjust.method,
                        group,
                        subject = NULL,
                        method = "ancombc2"){
  
  if(method == "ancombc2"){
    ## We add 0.5 as pseudo-count for ANCOM-BC2.
    assays <- S4Vectors::SimpleList(counts = as.matrix(feature.table + 0.5))
  }else if(method == "ancombc"){
    ## ANCOM-BC will add 1 as pseudo-count in it's own function.
    assays <- S4Vectors::SimpleList(counts = as.matrix(feature.table))
  }else{
    stop("Please use 'ancombc' or 'ancombc2' for analysis.")
  }
  
  smd <- S4Vectors::DataFrame(meta)
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
  if(!is.null(subject)){
    subject <- paste0("(1 | ", subject, " )")
  }
  if(method == "ancombc"){
    model = ancombc(data = tse, assay_name = "counts", 
                    tax_level = NULL, phyloseq = NULL, 
                    formula = formula,  
                    p_adj_method = adjust.method, 
                    prv_cut = 0, lib_cut = 0, 
                    group = group, struc_zero = FALSE, 
                    neg_lb = FALSE, tol = 1e-5, 
                    max_iter = 100, conserve = TRUE, 
                    alpha = 0.05, global = FALSE,
                    n_cl = 1, verbose = FALSE)
  }else if (method == "ancombc2"){
    model <- ancombc2(data = tse, assay_name = "counts", tax_level = NULL, 
                      fix_formula = formula, rand_formula = subject,
                      p_adj_method = adjust.method, pseudo = 0,
                      prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                      group = group, struc_zero = FALSE, neg_lb = FALSE,
                      alpha = 0.05, n_cl = 1, verbose = TRUE,
                      global = FALSE, pairwise = FALSE,
                      dunnet = FALSE, trend = FALSE)
                      # iter_control = list(tol = 1e-5, max_iter = 20,
                      #                     verbose = FALSE),
                      # em_control = list(tol = 1e-5, max_iter = 100))
  }
  return(model)
}