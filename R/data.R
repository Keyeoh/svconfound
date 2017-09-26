#' Small methylation example dataset
#'
#' A dataset containing a subset of a DNA Methylation analysis project based on
#' microarrays. It is used by a vignette in order to show how the package works.
#' Rows indicate probes and columns samples.
#'
#' @format A matrix with 3000 rows and 25 columns.
'smallmeth'

#' Partial phenotype information for smallmeth
#'
#' A data.frame describing part of the available phenotype information for the
#' methylation analysis project represented by the smallmeth dataset. Each row
#' describes a sample, and each column represents a different variable.
#'
#' @format A data.frame with 25 rows and 6 variables:
#' \describe{
#' \item{Sample_Group}{Group of the sample (Normal, Tumor)}
#' \item{Slide}{Slide in the plate where the sample was placed}
#' \item{Sex}{Gender of the sample (Male, Female)}
#' \item{DFS}{Disease Free Survival}
#' \item{GS}{Glutamine Synthetase}
#' \item{Age}{Age of the individual}
#' }
'phenopart'
