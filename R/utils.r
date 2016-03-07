#' Cut p-values according to a significance criterium.
#'
#' @param p_values A numeric containing the p-values
#'
#' @return A factor containing the cut categories.
cut_pvalues = function(p_values) {
  result = cut(p_values, c(0, 1e-10, 1e-5, 0.001, 0.05, Inf),
               labels = c('< 1e-10', '< 1e-5', '< 0.001', '< 0.05', 'NS'))
}

#' Compute associations with phenotypical variables.
#'
#' @param values A matrix containing either principal components or surrogate
#' variables.
#' @param pdata A data.frame containing the phenotypical data for the samples.
#' @param component_names A character vector containing the component names.
#'
#' @return A tbl_df with the results of the association tests.
#' @importFrom dplyr %>% mutate_
#' @importFrom tidyr gather_
compute_significance_data = function (values, pdata, component_names) {
  var_names = get_var_names(pdata)

  fits = lapply(var_names, function(xx) lm(values ~ pdata[[xx]]))
  sfits = lapply(fits, summary)
  pvalues = lapply(
    sfits,
    function(xx) sapply(xx, function(yy) get_p_value(yy$fstatistic))
    )
  pvalues_matrix = as.matrix(na.omit(t(as.data.frame(pvalues))))
  colnames(pvalues_matrix) = component_names

  significance_data = pvalues %>%
    as.data.frame() %>%
    mutate_('PC' = ~ component_names) %>%
    gather_('Variable', 'P_Value', var_names) %>%
    mutate_('Sig' = ~ cut_pvalues(P_Value))
}
