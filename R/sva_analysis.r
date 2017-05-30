#' Perform SVA confounder analysis
#'
#' Perform SVA confounder analysis and association tests between the Surrogate
#' Variables of the given input values and the phenotype data.
#'
#' This function performs a series of association tests between the Surrogate
#' Variables of the input values and the provided phenotype data.
#'
#' If a RGChannelSet is given, it can also compute the association with the
#' control probes. This implementation uses Linear models.
#'
#' @param values A matrix containing the data to analyze. Rows represent
#'   variables and columns samples.
#' @param pdata A data.frame containing the phenotype information for the
#'   samples.
#' @param main_formula A formula describing the main model of interest.
#' @param null_formula A formula describing the null model of interest.
#' @param vfilter The number of most variable rows to use in order to
#'   estimate the number of surrogate variables.
#' @param rgset If not NULL, compute association with control probes data
#'   contained in the provided RGChannelSet.
#' @return A list containing the number of Surrogate Variables (num_sv), the
#'   results of the confounding analysis tests (significance), and the Surrogate
#'   Variables obtained by SVA (surrogates).
#' @importFrom sva num.sv sva
#' @export
sva_analysis = function(values, pdata, main_formula, null_formula = NULL,
                        vfilter = NULL, rgset = NULL) {

  if (!is.null(vfilter)) {
    vfilter = min(vfilter, nrow(values))
  }

  mod = model.matrix(main_formula, data = pdata)
  mod0 = NULL

  if (!is.null(null_formula)) {
    mod0 = model.matrix(null_formula, data = pdata)
  }

  num_sv = num.sv(values, mod = mod, vfilter = vfilter)

  if (num_sv == 0) {
    stop(paste(
      'Number of estimated surrogate variables is zero.',
      'Check model and/or dataset.'
    ))
  }

  svs = sva(values, mod = mod, mod0 = mod0, n.sv = num_sv)

  sv_names = paste0('SV-', 1:(svs$n.sv))
  sv_names = factor(sv_names, levels = sv_names)

  if (num_sv == 1) {
    message(
      paste('If there is only a surrogate variable, the output from sva()',
            'is not ok, as it does not use the drop=FALSE parameter.',
            'Converting the surrogate variable back to a matrix.')
    )
    svs$sv = matrix(svs$sv, ncol = 1)
  }

  significance_data = compute_significance_data(svs$sv, pdata, sv_names)

  sig_data_rgset = NULL

  if (!is.null(rgset)) {
    data_control_values = get_control_variables(rgset)
    var_names = colnames(data_control_values)
    names(var_names) = var_names
    sig_data_rgset = compute_significance_data_var_names(
      svs$sv,
      data_control_values,
      sv_names,
      var_names
    )
    significance_data = rbind(significance_data, sig_data_rgset)

    significance_data$Variable = factor(
      significance_data$Variable,
      levels = unique(significance_data$Variable)
    )
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
