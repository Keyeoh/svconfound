#' Plot variance explained
#'
#' @param svd_result A list containing the result from the svd_analysis
#' function.
#'
#' @return A ggplot2 plot.
#' @importFrom ggplot2 ggplot aes_string geom_line geom_point xlab ylab theme_bw
#' theme element_text geom_vline
plot_svd_scree = function(svd_result) {
  result = ggplot(svd_result$variance_explained,
                  aes_string(x = 'PC', y = 'Var')) +
    geom_line(aes_string(group = '1'), col = 'red', linetype = 5) +
    geom_vline(xintercept = svd_result$limited_significant_PC, color = 'grey') +
    geom_point(size = 5, shape = 21, col = 'red', fill = 'lightblue') +
    ylab('Variance Explained') +
    xlab('Component') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  return(result)
}

#' Plot variance explained
#'
#' @param x A list containing the result from the svd_analysis
#' function.
#' @param ... Arguments to be passed to methods.
#' @return A ggplot2 plot.
#' @importFrom stats screeplot
#' @method screeplot SVDAnalysis
#' @export
screeplot.SVDAnalysis = function(x, ...) {
  return(plot_svd_scree(x, ...))
}

#' Plot results of confounder analysis.
#'
#' @param result A list containing the results from the confounder analysis
#' function.
#'
#' @return A ggplot2 plot.
#' @importFrom ggplot2 ggplot aes_string geom_tile scale_fill_manual
#' scale_color_manual xlab theme_bw theme geom_vline scale_x_discrete
plot_significance = function(result) {
  colors_palette = c('darkred', 'red', 'orange', 'gold', 'pink', 'transparent')
  result = ggplot(result$significance, aes_string(x = 'PC', y = 'Variable')) +
    geom_vline(xintercept = result$limited_significant_PC, color = 'grey') +
    scale_x_discrete() +
    geom_tile(aes_string(fill = 'Sig', color = 'Sig'), width = 0.6,
              height = 0.6) +
    scale_fill_manual(
      values = colors_palette,
      drop = FALSE,
      name = 'P-Value') +
    scale_color_manual(
      values = colors_palette,
      drop = FALSE,
      name = 'P-Value') +
    xlab('Component') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  return(result)
}

#' Plot results of confounder analysis.
#'
#' @param x A list containing the results from the confounder analysis
#' function.
#' @param ... Arguments to be passed to methods
#' @return A ggplot2 plot.
#' @method plot SVDAnalysis
#' @importFrom graphics plot
#' @export
plot.SVDAnalysis = function(x, ...) {
  return(plot_significance(x, ...))
}

#' Plot results of sv analysis.
#'
#' @param x A list containing the results from the sv analysis
#' function.
#' @param ... Arguments to be passed to methods
#' @return A ggplot2 plot.
#' @importFrom graphics plot
#' @method plot SVAAnalysis
#' @export
plot.SVAAnalysis = function(x, ...) {
  x$limited_significant_PC = x$num_sv
  return(plot_significance(x, ...))
}

