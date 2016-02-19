#' Perform SVA confounder analysis
#'
#' @param values A matrix containing the data to analyze. Rows are for variables and columns for samples.
#' @param pdata A data.frame containing the phenotypical information for the samples.
#' @param main_formula A formula describing the main model of interest.
#' @param vfilter The number of most variant variables to use in order to estimate the number of
#' surrogate variables.
#'
#' @return A list containing the number of surrogate variables (num_sv), the results of the confounding
#' analysis tests (significance), and the surrogate variables obtained by sva (surrogates).
#' @export
#' @importFrom sva num.sv sva
sva_analysis = function(values, pdata, main_formula, vfilter = 10000) {
  mod = model.matrix(main_formula, data = pdata)
  num_sv = num.sv(values, mod, vfilter = vfilter)
  svs = sva(values, mod, n.sv = num_sv)

  sv_names = paste0('SV-', 1:(svs$n.sv))
  sv_names = factor(sv_names, levels = sv_names)

  significance_data = compute_significance_data(svs$sv, pdata, sv_names)

  colnames(svs$sv) = gsub('-', '_', sv_names)

  result = list(
    num_sv = num_sv,
    significance = significance_data,
    surrogates = svs$sv
  )
}

#' Plot results of SVA confounder analysis.
#'
#' @param sva_result A list containing the results from the sva_analysis function.
#'
#' @return A ggplot2 plot.
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_tile scale_fill_manual scale_color_manual xlab theme_bw
#' theme
plot_sva_significance = function(sva_result) {
  result = ggplot(sva_result$significance, aes_string(x = 'PC', y = 'Variable')) +
    geom_tile(aes_string(fill = 'Sig', color = 'Sig'), width = 0.6, height = 0.6) +
    scale_fill_manual(values = c('darkred', 'red', 'orange', 'pink', 'transparent'), drop = FALSE,
                      name = 'P-Value') +
    scale_color_manual(values = c('darkred', 'red', 'orange', 'pink', 'transparent'), drop = FALSE,
                       name = 'P-Value') +
    xlab('Component') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  return(result)
}
