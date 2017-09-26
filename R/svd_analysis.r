#' Perform SVD confounder analysis.
#'
#' Perform SVD confounder analysis and association tests between the Principal
#' Components of the given input values and the phenotype data.
#'
#' This function performs a series of association tests between the Principal
#' Components of the input values and the provided phenotype data. Previously,
#' it centers and scales the data, something that might be useful when the
#' variables are measured in different units.
#'
#' If a RGChannelSet is given, it can also compute the association with the
#' control probes. This implementation can use two different methods for testing
#' (Linear models and Kruskal-Wallis).
#'
#' The number of meaningful components is computed using the function
#' \code{EstDimRMT} from the package \code{isva}.
#'
#' @param values A matrix of numerical values. Rows represent
#'   variables and columns samples.
#' @param pdata A data.frame containing the phenotype data for the samples.
#' @param center A logical value indicating whether the variables should be
#'   shifted to be zero centered. Alternately, a vector of length equal the
#'   number of columns of x can be supplied. The value is passed to scale.
#' @param scale A logical value indicating whether the variables should be
#'   scaled to have unit variance before the analysis takes place. The default
#'   is FALSE for consistency with S, but in general scaling is advisable.
#'   Alternatively, a vector of length equal the number of columns of x can be
#'   supplied. The value is passed to scale. Use when the variables are in
#'   arbitrary units of measurement.
#' @param rgset If not NULL, compute association with control probes data
#'   contained in the provided RGChannelSet.
#' @param method Method used for computing p-values.
#' @return A list containing the proportion of variance explained by the
#'   principal components (variance_explained) and a data.frame representing the
#'   results from the association analysis (significance).
#'
#' @export
svd_analysis = function(
  values,
  pdata,
  center = TRUE,
  scale = FALSE,
  rgset = NULL,
  method = c('lm', 'kruskal')) {

  method = match.arg(method)
  scaled_values = scale(t(values), center = center, scale = scale)

  sv_decomp = svd(scaled_values)
  var_explained = sv_decomp[['d']] ^ 2 / sum(sv_decomp[['d']] ^ 2)
  component_names = paste0('PC-', 1:length(var_explained))
  component_names = factor(component_names, levels = component_names)

  variance_explained_data = data.frame(
    PC = component_names,
    Var = var_explained
  )

  if (nrow(pdata) != nrow(scaled_values)) {
    stop(paste('The num of pdata row elements must be equal',
               'to row elements of [v] or [u].',
               'Actual number of row elements pdata:', nrow(pdata)))
  }

  significance_data = compute_significance_data(sv_decomp[['u']], pdata,
                                                component_names,
                                                method = method)

  sig_data_rgset = NULL

  if (!is.null(rgset)) {
    data_control_values = get_control_variables(rgset)
    var_names = colnames(data_control_values)
    names(var_names) = var_names
    sig_data_rgset = compute_significance_data_var_names(sv_decomp[['u']],
                                                         data_control_values,
                                                         component_names,
                                                         var_names)
    significance_data = rbind(significance_data, sig_data_rgset)

    significance_data[['Variable']] = factor(
      significance_data[['Variable']],
      levels = unique(significance_data[['Variable']])
      )
  }

  dim_pca = isva::EstDimRMT(t(scaled_values), plot = FALSE)[['dim']]

  result = list(
    variance_explained = variance_explained_data,
    center = center,
    scale = scale,
    method = method,
    significance = significance_data,
    significance_control = sig_data_rgset,
    limit_significant_PC = dim_pca
  )

  class(result) = append(class(result), 'SVDAnalysis')
  return(result)
}
