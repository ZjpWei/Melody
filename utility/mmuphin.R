fit_metaAnalysis <- function(feature.abd,
                             data,
                             test_variable,
                             contrasts,
                             batch_variable,
                             covariates = NULL,
                             covariates.random = NULL,
                             moderator_variables = NULL,
                             n_features = 50,
                             normalization = "TSS",
                             transform = "AST",
                             analysis_method = "LM",
                             rma.method = "REML",
                             forest.plots = TRUE,
                             verbose = TRUE,
                             rma.threshold = 1e-6,
                             rma.maxiter = 1000) {
 
  result <- MMUPHin::lm_meta(
    feature_abd = feature.abd,
    batch = batch_variable,
    exposure = test_variable,
    covariates = covariates,
    covariates_random = covariates.random,
    data = data,
    control = list(
      normalization = normalization,
      transform = transform,
      analysis_method = analysis_method,
      rma_method = rma.method,
      verbose = verbose,
      rma_conv = rma.threshold,
      rma_maxit = rma.maxiter
    )
  )
  if(!is.null(moderator_variables)) {
    data.moderator <- data %>%
      dplyr::select(dplyr::one_of(c(batch_variable, moderator_variables))) %>%
      dplyr::group_by(!!sym(batch_variable)) %>%
      dplyr::summarise_all(function(x) unique(x)) %>%
      dplyr::mutate_at(moderator_variables, factor) %>%
      as.data.frame() %>%
      tibble::column_to_rownames(batch_variable)
    df_moderator <- MMUPHin:::rma_wrapper(result$maaslin_fits,
                                          # data.moderator = data.moderator,
                                          method = rma.method,
                                          # rma.threshold = rma.threshold,
                                          rma_maxit = rma.maxiter)
  }
  
  return(result)
}

