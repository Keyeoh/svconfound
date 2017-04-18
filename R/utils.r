#' Cut p-values.
#'
#' Cut p-values according to a significance criterium.
#'
#' This is a simple internal function for cutting a vector of p-values into
#' a set of meaningful categories. The categories are taken from the original
#' book chapter that inspired the package.
#'
#' @param p_values A numeric containing the p-values
#' @return A factor containing the resulting categories.
cut_pvalues = function(p_values) {
  result = cut(
    p_values, c(-Inf, 1e-10, 1e-5, 0.001, 0.01, 0.05, Inf),
    labels = c('< 1e-10', '< 1e-5', '< 0.001', '< 0.01', '< 0.05', 'NS')
    )
  return(result)
}

#' Get p-value from f statistic data.
#'
#' @param fstat A numeric vector of length three obtained from a linear model
#' fit.
#'
#' @return The computed p-value.
get_p_value = function(fstat) {
  if (is.numeric(fstat['value'])) {
    p_value = 1 - pf(fstat['value'], fstat['numdf'], fstat['dendf'])
  } else {
    p_value = 1
  }
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

#' Compute associations with phenotypical variables.
#'
#' Compute associations with phenotypical and control variables in order to
#' find out which of them are significant, and thus possibly introducing
#' some confounding in the data.
#'
#' The current implementation allows the user to use two different methods
#' for testing association: a general linear model (lm) and a Kruskal-Wallis
#' test. User should be careful, as the Kruskal Wallis test should not be used
#' for continuous independent variables.
#'
#' An association test is performed between every component (Principal Component
#' or Surogate Variable) and every phenotype and control variable. P-values are
#' then collected and stored in a \code{tbl_df} object describing all the
#' comparisons.
#'
#' This function receives variable names as a parameter.
#'
#' @inheritParams compute_significance_data
#' @param var_names Variable names to be used for calculating the p_values.
#' @return A tbl_df with the results of the association tests.
#' @importFrom dplyr %>% mutate_
#' @importFrom tidyr gather_
#' @importFrom stats lm na.omit pf model.matrix kruskal.test
compute_significance_data_var_names = function(values, pdata, component_names,
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

  return(significance_data)
}

#' Compute associations with phenotypical variables.
#'
#' Compute associations with phenotypical and control variables in order to
#' find out which of them are significant, and thus possibly introducing
#' some confounding in the data.
#'
#' The current implementation allows the user to use two different methods
#' for testing association: a general linear model (lm) and a Kruskal-Wallis
#' test. User should be careful, as the Kruskal Wallis test should not be used
#' for continuous independent variables.
#'
#' An association test is performed between every component (Principal Component
#' or Surogate Variable) and every phenotype and control variable. P-values are
#' then collected and stored in a \code{tbl_df} object describing all the
#' comparisons.
#'
#' @param values A matrix containing either principal components or surogate
#' variables.
#' @param pdata A data.frame containing the phenotypical data for the samples.
#' @param component_names A character vector containing the component names.
#' @param method A character indicating the method used for the association
#' tests.
#' @return A tbl_df with the results of the association tests.
#' @importFrom dplyr %>% mutate_
#' @importFrom tidyr gather_
compute_significance_data = function (values, pdata, component_names,
                                      method = c('lm', 'kruskal')) {
  var_names = get_var_names(pdata)

  return(compute_significance_data_var_names(values, pdata, component_names,
                                        var_names, method))
}

#' Get control variables.
#'
#' Get a matrix containing the signals for the control probes of a given
#' RGChannelSet object.
#'
#' The function \code{get_control_variables} just takes a \code{RGChannelSet}
#' object as input, accesses the values of its control probes, and returns
#' them as a numeric matrix.
#'
#' This is an internal function. Its purpose is to provide the main analysis
#' functions with more data to test. It is usually useful to test if control
#' probes are introducing some confounding in the dataset, as it might be
#' signaling that there are unknown effects in the data.
#'
#' @param rgset An RGChannelSet object containing the control data.
#' @importFrom methods is
#' @importFrom minfi getControlAddress getRed getGreen
#' @return A matrix with the values of the control elements.
get_control_variables = function(rgset) {
  datac2_m = NULL

  if (!is.null(rgset) && is(rgset, 'RGChannelSet')) {
    bc1 = getControlAddress(rgset, controlType = c('BISULFITE CONVERSION I'))
    bc2 = getControlAddress(rgset, controlType = c('BISULFITE CONVERSION II'))
    ext = getControlAddress(rgset, controlType = c('EXTENSION'))
    tr = getControlAddress(rgset, controlType = c('TARGET REMOVAL'))
    hyb = getControlAddress(rgset, controlType = c('HYBRIDIZATION'))

    control_names = c(
      'BSC-I C1 Grn',
      'BSC-I C2 Grn',
      'BSC-I C3 Grn',
      'BSC-I C4 Red',
      'BSC-I C5 Red',
      'BSC-I C6 Red',
      'BSC-II C1 Red',
      'BSC-II C2 Red',
      'BSC-II C3 Red',
      'BSC-II C4 Red',
      'Target Removal 1 Grn',
      'Target Removal 2 Grn',
      'Hyb (Low) Grn',
      'Hyb (Medium) Grn',
      'Hyb (High) Grn',
      'Extension (A) Red',
      'Extension (T) Red',
      'Extension (C) Grn',
      'Extension (G) Grn'
    )

    control_signals = rbind(
      getGreen(rgset)[bc1[1:3],],
      getRed(rgset)[bc1[7:9],],
      getRed(rgset)[bc2[1:4],],
      getGreen(rgset)[tr[1:2],],
      getGreen(rgset)[hyb[1:3],],
      getRed(rgset)[ext[1:2],],
      getGreen(rgset)[ext[3:4],]
    )
    dimnames(control_signals) = list(control_names, colnames(rgset))

    datac2_m = t(log2(control_signals))
  }

  return(datac2_m)
}
