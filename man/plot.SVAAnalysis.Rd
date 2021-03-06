% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.r
\name{plot.SVAAnalysis}
\alias{plot.SVAAnalysis}
\title{Plot results of SVA confounder analysis.}
\usage{
\method{plot}{SVAAnalysis}(x, ...)
}
\arguments{
\item{x}{A list containing the results from the confounder analysis
function.}

\item{...}{Arguments to be passed to methods}
}
\value{
A ggplot2 heatmap showing the results of all the association tests.
}
\description{
Plot results of SVA confounder analysis as a heatmap where each cell
represents the statistical significance of the given association test.
}
\details{
This internal function takes the result of the \code{sva_analysis} function
and returns a heatmap representing the degree of significance for each of the
performed association tests. Surrogate Variables are represented on the x
axis. Control and phenotype variables are represented in the y axis. The
color of each cell corresponds to the degree of statistical significance for
the corresponding association test.
}
