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

#' Plot variance explained
#'
#' @param svd_result A list containing the result from the svd_analysis function.
#'
#' @return A ggplot2 plot.
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_line geom_point xlab ylab theme_bw theme element_text
plot_svd_scree = function(svd_result) {
  result = ggplot(svd_result$variance_explained, aes_string(x = 'PC', y = 'Var')) +
    geom_line(aes_string(group = '1'), col = 'red', linetype = 5) +
    geom_point(size = 5, shape = 21, col = 'red', fill = 'lightblue') +
    ylab('Variance Explained') +
    xlab('Component') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  return(result)
}

#' Plot results of SVD confounder analysis.
#'
#' @param svd_result A list containing the result from the svd_analysis function.
#'
#' @return A ggplot2 plot.
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_tile scale_fill_manual xlab theme_bw theme
plot_svd_significance = function(svd_result) {
  result = ggplot(svd_result$significance, aes_string(x = 'PC', y = 'Variable')) +
    geom_tile(aes_string(fill = 'Sig')) +
    scale_fill_manual(values = c('darkred', 'red', 'orange', 'pink', 'white'), drop = FALSE,
                      name = 'P-Value') +
    xlab('Component') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  return(result)
}
