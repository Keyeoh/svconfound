#' Cut p-values according to a significance criterium.
#'
#' @param p_values A numeric containing the p-values
#'
#' @return A factor containing the cut categories.
cut_pvalues = function(p_values) {
  result = cut(p_values, c(0, 1e-10, 1e-5, 0.001, 0.05, Inf),
               labels = c('< 1e-10', '< 1e-5', '< 0.001', '< 0.05', 'NS'))
}
