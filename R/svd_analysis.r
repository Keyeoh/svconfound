#' Get p-value from f statistic data.
#'
#' @param fstat A numeric vector of length three obtained from a linear model
#' fit.
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
  var_names = colnames(pdata)[(n_levels > 1 & n_levels < nrow(pdata)) |
                                sapply(pdata, is.numeric)]
  names(var_names) = var_names

  return(var_names)
}

#' Get control variables
#'
#' @param rgSet. RGSet data with control elements
#'
#' @return A matrix with the values of the control elements.
get_control_variables = function(rgSet) {
  dataC2.m = NULL
  if (!is.null(rgSet)) {
    bc1 = getControlAddress(rgSet, controlType = c('BISULFITE CONVERSION I'))
    bc2 = getControlAddress(rgSet, controlType = c('BISULFITE CONVERSION II'))
    ext = getControlAddress(rgSet, controlType = c('EXTENSION'))
    tr = getControlAddress(rgSet, controlType = c('TARGET REMOVAL'))
    hyb = getControlAddress(rgSet, controlType = c('HYBRIDIZATION'))
    CPP = rbind(getGreen(rgSet)[bc1[1:3], ],
                getRed(rgSet)[bc1[7:9], ],
                getRed(rgSet)[bc2[1:4], ],
                getGreen(rgSet)[tr[1:2], ],
                getGreen(rgSet)[hyb[1:3], ],
                getRed(rgSet)[ext[1:2], ],
                getGreen(rgSet)[ext[3:4], ])
    controlNames = c('BSC-I C1 Grn', 'BSC-I C2 Grn', 'BSC-I C3 Grn',
                      'BSC-I C4 Red', 'BSC-I C5 Red', 'BSC-I C6 Red',
                      'BSC-II C1 Red', 'BSC-II C2 Red', 'BSC-II C3 Red',
                      'BSC-II C4 Red', 'Target Removal 1 Grn', 'Target Removal 2 Grn',
                      'Hyb (Low) Grn', 'Hyb (Medium) Grn', 'Hyb (High) Grn',
                      'Extension (A) Red', 'Extension (T) Red', 'Extension (C) Grn',
                      'Extension (G) Grn')
    rownames(CPP) = controlNames
    colnames(CPP) = colnames(rgSet)
    dataC2.m = t(log2(CPP))
    # matrix(nrow = length(controlNames), ncol = ncol(rgSet))
    rownames(dataC2.m) = colnames(rgSet)
    colnames(dataC2.m) = controlNames
#     for (r in 1:nrow(dataC2.m)) {
#       dataC2.m[r, ] = log2(as.numeric(CPP[r, ]))
#     }
  }
  return(dataC2.m)
}

#' Compute SVD confounder analysis.
#'
#' @param values A matrix of numerical values. Rows represent variables and
#' columns samples.
#' @param pdata A data.frame containing the phenotypical data for the samples.
#' #' @param center A logical value indicating whether the variables should be shifted to
#' be zero centered. Alternately, a vector of length equal the number of columns of x
#' can be supplied. The value is passed to scale.
#' @param scale A logical value indicating whether the variables should be scaled to
#'  have unit variance before the analysis takes place. The default is FALSE for
#'   consistency with S, but in general scaling is advisable. Alternatively, a
#'   vector of length equal the number of columns of x can be supplied. The value is
#'   passed to scale. Use when the variables are in arbitrary units of measurement
#' @param significante_data_from_samples. If TRUE obtain p-values and variance explained from samples.
#' Otherwise obtain it from variables. Change only experts for obtaining other data
#' @param rgSet If not NULL, calculate show control variables
#' @param method, Method used for calculate p-values for pdata variables
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
               'Actual number of row elements pdata:', num_samples))
  }

  significance_data = compute_significance_data(matrix_for_sig_data, pdata,
                                                component_names,
                                                method = method)

  data_control_values = get_control_variables(rgSet)

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
    dim_pca = isva::EstDimRMT(t(values[1:max_samples, ]))$dim
  } else {
    # Rows label features/variables, Columns samples.
    dim_pca = isva::EstDimRMT(t(values))$dim
  }

  result = list(
    variance_explained = variance_explained_data,
    significance = significance_data,
    significance_control = sig_data_rgset,
    limited_significant_PC = dim_pca
  )

  class(result) = append(class(result), 'SVDAnalysis')

  # return(structure(result, class = 'SVDAnalysis'))
  return(result)
}
