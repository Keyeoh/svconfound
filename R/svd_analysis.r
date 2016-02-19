#' Get p-value from f statistic data.
#'
#' @param fstat A numeric vector of length three obtained from a linear model fit.
#'
#' @return The computed p-value.
get_p_value = function(fstat) {
  p_value = 1 - pf(fstat['value'], fstat['numdf'], fstat['dendf'])
  return(p_value)
}

#' Get variable names suitable for analysis.
#'
#' @param pdata A data.frame containing the phenotypical data.
#'
#' @return A character vector containing the selected variable names.
get_var_names = function(pdata) {
  n_levels = sapply(pdata, function(xx) length(unique(xx)))
  var_names = colnames(pdata)[n_levels > 1 & n_levels < 9]
  names(var_names) = var_names

  return(var_names)
}

#' Compute SVD confounder analysis.
#'
#' @param values A matrix of numerical values. Rows represent variables and columns samples.
#' @param pdata A data.frame containing the phenotypical data for the samples.
#'
#' @return A list containing the proportion of variance explained by the principal components
#' (variance_explained) and a data.frame representing the results from the association analysis
#' (significance).
#'
#' @export
svd_analysis = function(values, pdata) {
  sv_decomp = svd(values)
  var_explained = sv_decomp$d ^ 2 / sum(sv_decomp$d ^ 2)
  component_names = paste0('PC-', 1:length(var_explained))
  component_names = factor(component_names, levels = component_names)

  variance_explained_data = data.frame(
    PC = component_names,
    Var = var_explained
  )

  significance_data = compute_significance_data(sv_decomp$v, pdata, component_names)

  result = list(
    variance_explained = variance_explained_data,
    significance = significance_data
  )

  return(result)
}
