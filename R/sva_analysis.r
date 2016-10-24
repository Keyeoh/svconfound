#' Perform SVA confounder analysis
#'
#' @param values A matrix containing the data to analyze. Rows are for variables
#' and columns for samples.
#' @param pdata A data.frame containing the phenotypical information for the
#' samples.
#' @param main_formula A formula describing the main model of interest.
#' @param null_formula A formula describing the null model of interest.
#' @param vfilter The number of most variant variables to use in order to
#' estimate the number of surrogate variables.
#'
#' @return A list containing the number of surrogate variables (num_sv), the
#' results of the confounding analysis tests (significance), and the surrogate
#' variables obtained by sva (surrogates).
#'
#' @export
#' @importFrom sva num.sv sva
sva_analysis = function(values, pdata, main_formula,
                        null_formula = NULL,
                        vfilter = NULL,
                        rgSet = NULL) {

  if (!is.null(vfilter)) {
    vfilter = min(vfilter, nrow(values))
  }

  mod = model.matrix(main_formula, data = pdata)
  mod0 = NULL
  if (!is.null(null_formula)) {
    mod0 = model.matrix(null_formula, data = pdata)
  }
  num_sv = num.sv(values, mod = mod, vfilter = vfilter)
  svs = sva(values, mod = mod, mod0 = mod0, n.sv = num_sv)

  sv_names = paste0('SV-', 1:(svs$n.sv))
  sv_names = factor(sv_names, levels = sv_names)

  significance_data = compute_significance_data(svs$sv, pdata, sv_names)

  data_control_values = get_control_variables(rgSet)
  sig_data_rgset = NULL

  if (!is.null(data_control_values)) {
    var_names = colnames(data_control_values)
    names(var_names) = var_names
    sig_data_rgset = compute_significance_data_var_names(svs$sv,
                                                         data_control_values,
                                                         sv_names,
                                                         var_names)
    significance_data = rbind(significance_data, sig_data_rgset)
    # for sorting elements
    significance_data$Variable = factor(significance_data$Variable,
                                        levels = unique(significance_data$Variable))
  }

  colnames(svs$sv) = gsub('-', '_', sv_names)

  result = list(
    num_sv = num_sv,
    mod0 = mod0,
    mod = mod,
    significance = significance_data,
    significance_control = sig_data_rgset,
    surrogates = svs$sv
  )

  class(result) = append(class(result), 'SVAAnalysis')
  return(result)
}
