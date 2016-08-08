#' Cut p-values according to a significance criterium.
#'
#' @param p_values A numeric containing the p-values
#'
#' @return A factor containing the cut categories.
cut_pvalues = function(p_values) {
  result = cut(p_values, c(0, 1e-10, 1e-5, 0.001, 0.01, 0.05, Inf),
               labels = c('< 1e-10', '< 1e-5', '< 0.001', '< 0.01', '< 0.05',
                          'NS'))
}


#' Compute associations with phenotypical variables.
#'
#' @param values A matrix containing either principal components or surrogate
#' variables.
#' @param pdata A data.frame containing the phenotypical data for the samples.
#' @param component_names A character vector containing the component names.
#' @param var_names Variable names to be used for calculatin de p_values
#' @param method. The algoritm used for calculating p-values
#'
#' @return A tbl_df with the results of the association tests.
#' @importFrom dplyr %>% mutate_
#' @importFrom tidyr gather_
compute_significance_data_var_names = function (values, pdata, component_names,
                                                var_names,
                                                method = c('lm', 'kruskal')) {
  method = match.arg(method)
  if (method == 'kruskal') {
    fits = function(val, pd) {
      kruskal.test(val, pd)$p.value
    }
    pvalues = lapply(var_names,
                     function(xx) apply(values, 2,
                                        fits,
                                        as.factor(pdata[, xx])))
  } else {
    fits = lapply(var_names, function(xx) lm(values ~ pdata[, xx]))
    sfits = lapply(fits, summary)
    pvalues = lapply(
      sfits,
      function(xx) sapply(xx, function(yy) get_p_value(yy$fstatistic))
    )
  }
  p_values_df = data.frame(pvalues, check.names = FALSE)

  significance_data = p_values_df %>%
    mutate_('PC' = ~ component_names) %>%
    gather_('Variable', 'P_Value', var_names) %>%
    mutate_('Sig' = ~ cut_pvalues(P_Value))
}

#' Compute associations with phenotypical variables.
#'
#' @param values A matrix containing either principal components or surrogate
#' variables.
#' @param pdata A data.frame containing the phenotypical data for the samples.
#' @param component_names A character vector containing the component names.
#' @param method. The algoritm used for calculating p-values
#'
#' @return A tbl_df with the results of the association tests.
#' @importFrom dplyr %>% mutate_
#' @importFrom tidyr gather_
compute_significance_data = function (values, pdata, component_names,
                                      method = c('lm', 'kruskal')) {
  var_names = get_var_names(pdata)
  return(compute_significance_data_var_names(values, pdata, component_names, var_names,
                                             method))
}
