#' Plot ratio of explained variance.
#'
#' Generates a plot showing the ratio of explained variance by the Principal
#' Components as a scree plot.
#'
#' This is an internal function taking the output of the \code{svd_analysis}
#' function as input. It generates a scree plot showing the ratio of explained
#' variance for the ordered Principal Components. It is very similar to the
#' \code{screeplot} function but it uses ggplot2 instead of base graphics.
#'
#' @param x A list containing the results from the svd_analysis
#' function.
#' @return A ggplot2 scree plot showing the ratio of explained variance for
#' every Principal Component.
#' @importFrom ggplot2 ggplot aes_string geom_line geom_point xlab ylab theme_bw
#' theme element_text geom_vline
#' @keywords internal
plot_svd_scree = function(x) {
  result = ggplot(x[['variance_explained']],
                  aes_string(x = 'PC', y = 'Var')) +
    geom_line(aes_string(group = '1'), col = 'red', linetype = 5) +
    geom_vline(xintercept = x[['limit_significant_PC']], color = 'grey') +
    geom_point(size = 5, shape = 21, col = 'red', fill = 'lightblue') +
    ylab('Variance Explained') +
    xlab('Component') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  return(result)
}

#' Plot ratio of explained variance.
#'
#' Generates a plot showing the ratio of explained variance by the Principal
#' Components as a scree plot.
#'
#' This is a generic method taking the output of the \code{svd_analysis}
#' function as input. It generates a scree plot showing the ratio of explained
#' variance for the ordered Principal Components. It is very similar to the
#' \code{screeplot} function but it uses ggplot2 instead of base graphics.
#'
#' @param x A list containing the results from the svd_analysis function.
#' @param ... Arguments to be passed to methods.
#' @return A ggplot2 scree plot showing the ratio of explained variance for
#' every Principal Component.
#' @importFrom stats screeplot
#' @method screeplot SVDAnalysis
#' @export
screeplot.SVDAnalysis = function(x, ...) {
  return(plot_svd_scree(x, ...))
}

#' Plot ratio of explained variance.
#'
#' Generates an error when the user tries to call the \code{screeplot} function
#' on a SVAAnalysis object.
#'
#' This is a generic method taking the output of the \code{sva_analysis}
#' function as input. It generates an error..
#'
#' @param x A list containing the results from the sva_analysis function.
#' @param ... Arguments to be passed to methods.
#' @importFrom stats screeplot
#' @method screeplot SVAAnalysis
#' @export
screeplot.SVAAnalysis = function(x, ...) {
  stop('cannot be called on a SVAAnalysis object')
}

#' Plot results of SVD/SVA confounder analysis.
#'
#' Plot results of SVD/SVA confounder analysis as a heatmap where each cell
#' represents the statistical significance of the given association test.
#'
#' This internal function takes the result of either \code{svd_analysis} or
#' \code{sva_analysis} functions and returns a heatmap representing the degree
#' of significance for each of the performed association tests. Components or
#' surrogate variables are represented on the x axis. Control and phenotype
#' variables are represented in the y axis. The color of each cell corresponds
#' to the degree of statistical significance for the corresponding association
#' test.
#'
#' @param x A list containing the results from the confounder analysis
#' function.
#' @return A ggplot2 heatmap showing the results of all the association tests.
#' @importFrom ggplot2 ggplot aes_string geom_tile scale_fill_manual
#' scale_color_manual xlab theme_bw theme geom_vline scale_x_discrete
#' @keywords internal
plot_significance = function(x) {
  colors_palette = c('darkred', 'red', 'orange', 'gold', 'pink', 'transparent')

  x = ggplot(x[['significance']], aes_string(x = 'PC', y = 'Variable')) +
    geom_vline(xintercept = x[['limit_significant_PC']], color = 'grey') +
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

  return(x)
}

#' Plot results of SVD confounder analysis.
#'
#' Plot results of SVD confounder analysis as a heatmap where each cell
#' represents the statistical significance of the given association test.
#'
#' This internal function takes the result of the \code{svd_analysis} function
#' and returns a heatmap representing the degree of significance for each of the
#' performed association tests. Principal Components are represented on the x
#' axis. Control and phenotype variables are represented in the y axis. The
#' color of each cell corresponds to the degree of statistical significance for
#' the corresponding association test.
#'
#' @param x A list containing the results from the confounder analysis
#' function.
#' @param ... Arguments to be passed to methods
#' @return A ggplot2 heatmap showing the results of all the association tests.
#' @method plot SVDAnalysis
#' @importFrom graphics plot
#' @export
plot.SVDAnalysis = function(x, ...) {
  return(plot_significance(x, ...))
}

#' Plot results of SVA confounder analysis.
#'
#' Plot results of SVA confounder analysis as a heatmap where each cell
#' represents the statistical significance of the given association test.
#'
#' This internal function takes the result of the \code{sva_analysis} function
#' and returns a heatmap representing the degree of significance for each of the
#' performed association tests. Surrogate Variables are represented on the x
#' axis. Control and phenotype variables are represented in the y axis. The
#' color of each cell corresponds to the degree of statistical significance for
#' the corresponding association test.
#'
#' @param x A list containing the results from the confounder analysis
#' function.
#' @param ... Arguments to be passed to methods
#' @return A ggplot2 heatmap showing the results of all the association tests.
#' @method plot SVAAnalysis
#' @importFrom graphics plot
#' @export
plot.SVAAnalysis = function(x, ...) {
  x[['limit_significant_PC']] = x[['num_sv']]
  return(plot_significance(x, ...))
}

