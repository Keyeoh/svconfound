#' Compute SVD confounder analysis.
#'
#' @param values A matrix of numerical values. Rows represent samples and
#' columns variables.
#' @param pdata A data.frame containing the phenotypical data for the samples.
#' @param center A logical value indicating whether the variables should be shifted to
#' be zero centered. Alternately, a vector of length equal the number of columns of x
#' can be supplied. The value is passed to scale.
#' @param scale A logical value indicating whether the variables should be scaled to
#'  have unit variance before the analysis takes place. The default is FALSE for
#'   consistency with S, but in general scaling is advisable. Alternatively, a
#'   vector of length equal the number of columns of x can be supplied. The value is
#'   passed to scale. Use when the variables are in arbitrary units of measurement
#' @param significante_data_from_samples If TRUE obtain p-values and variance explained from samples.
#' Otherwise obtain it from variables. Change only experts for obtaining other data
#' @param rgSet If not NULL, calculate show control variables
#' @param method Method used for calculate p-values for pdata variables
#'
#' @return A list containing the proportion of variance explained by the
#' principal components (variance_explained) and a data.frame representing the
#' results from the association analysis (significance).
#'
#' @export
svd_analysis = function(values, pdata, center = T, scale = F,
                        significante_data_from_samples = nrow(pdata) == nrow(values),
                        rgSet = NULL,
                        method = c('lm', 'kruskal')) {

  # Rows label samples, Columns features/variables.
  values = scale(values, center = center, scale = scale)

  sv_decomp = svd(values)
  var_explained = sv_decomp$d ^ 2 / sum(sv_decomp$d ^ 2)
  component_names = paste0('PC-', 1:length(var_explained))
  component_names = factor(component_names, levels = component_names)

  variance_explained_data = data.frame(
    PC = component_names,
    Var = var_explained
  )

  # check the number of elements you want to test
  if (significante_data_from_samples) {
    matrix_for_sig_data = sv_decomp$u
  } else if (nrow(pdata) == ncol(values)) {
    matrix_for_sig_data = sv_decomp$v
  } else {
    stop(paste('The num of pdata row elements must be equal to row elements of [v] or [u].',
               'Actual number of row elements pdata:', nrow(pdata)))
  }

  significance_data = compute_significance_data(matrix_for_sig_data, pdata,
                                                component_names,
                                                method = method)

  data_control_values = get_control_variables(rgSet)
  sig_data_rgset = NULL

  if (!is.null(data_control_values)) {
    var_names = colnames(data_control_values)
    names(var_names) = var_names
    sig_data_rgset = compute_significance_data_var_names(matrix_for_sig_data,
                                                         data_control_values,
                                                         component_names,
                                                         var_names)
    significance_data = rbind(significance_data, sig_data_rgset)
    # for sorting elements
    significance_data$Variable = factor(significance_data$Variable,
                                        levels = unique(significance_data$Variable))
  }

  max_samples = 5000
  # the matrix row values indicates the samples and the columns the variables
  if (nrow(values) > max_samples) {
    # select the 5000 first samples
    # warning!! we are removing samples (the original values matrix has samples like rows
    # and variables like colums)
    # Rows label features/variables, Columns samples.
    dim_pca = isva::EstDimRMT(t(values[1:max_samples, ]), plot = FALSE)$dim
  } else {
    # Rows label features/variables, Columns samples.
    dim_pca = isva::EstDimRMT(t(values), plot = FALSE)$dim
  }

  result = list(
    variance_explained = variance_explained_data,
    center = center,
    scale = scale,
    method = method,
    significance = significance_data,
    significance_control = sig_data_rgset,
    limited_significant_PC = dim_pca
  )

  class(result) = append(class(result), 'SVDAnalysis')
  return(result)
}
