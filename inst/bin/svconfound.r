#!/usr/bin/env Rscript

library(docopt)
library(ggplot2)

# -------------------------------------------------------------------------------------------------
script_version = '0.1.0\n'

# -------------------------------------------------------------------------------------------------
doc = '
Usage:
  svconfound.r [--csv] [--prefix=PREFIX] [--fixed=VAR]... [--model=MODEL] <input_file>
  svconfound.r -h
  svconfound.r --version

Options:
  -c, --csv  Input file is in csv format.
  -p PREFIX, --prefix=PREFIX  Set output report identifier to NAME.
  -f VAR, --fixed=VAR  Sets a variable as fixed by the experimentalist.
  -m MODEL, --model=MODEL  Sets the model to be used by the SVA step.
  -h, --help  Show this message.
  --version  Print the version.
'

# -------------------------------------------------------------------------------------------------
perform_analysis = function(values, pdata, model, fixed_vars = NULL) {
  library(sva)
  library(svconfound)

  results = list()

  # Step 1. SVD analysis.
  results$svd_results = svd_analysis(values, pdata)

  # Step 2. ComBat (optional).
  if (!is.null(fixed_vars)) {
    batch_variable = Reduce(function(xx, yy) paste(xx, yy, sep = '_'),
                            pdata[, fixed_vars, drop = FALSE])
    fixed_values = ComBat(values, batch_variable)
  } else {
    fixed_values = values
  }

  # Step 3. SVA analysis.
  results$sva_results = sva_analysis(fixed_values, pdata, as.formula(model))

  return(results)
}

# -------------------------------------------------------------------------------------------------
if(!interactive()) {
  args = docopt(doc, version = script_version)

  if (is.null(args$prefix)) {
    args$prefix = strsplit(basename(args$input_file), split = '\\.')[[1]][1]
  }
  message(paste('Setting prefix to', args$prefix))

  if (is.null(args$model)) {
    args$model = '~ 1'
  }
  message(paste('Setting model to', args$model))

  if (args$csv) {
    stop('Not implemented yet!')
  } else {
    rds_object = tryCatch(readRDS(args$input_file),
                          error = function(e) {
                            stop('Input file is not a valid RDS file.')
                          })

    if (!inherits(rds_object, c('MethylSet', 'GenomicMethylSet', 'RatioSet', 'GenomicRatioSet'))) {
      stop('Input file must contain a minfi methylation object.')
    }

    values = getM(rds_object)
    pheno_data = pData(rds_object)
  }

  analysis_results = perform_analysis(values, pheno_data, args$model, args$fixed)

  write.table(analysis_results$svd_results$variance_explained,
              file = paste0(args$prefix, '.variance_explained.tsv'),
              sep = '\t', dec = '.', quote = FALSE, col.names = TRUE, row.names = FALSE)

  write.table(analysis_results$svd_results$significance,
              file = paste0(args$prefix, '.svd_significance.tsv'),
              sep = '\t', dec = '.', quote = FALSE, col.names = TRUE, row.names = FALSE)

  write.table(analysis_results$sva_results$significance,
              file = paste0(args$prefix, '.sva_significance.tsv'),
              sep = '\t', dec = '.', quote = FALSE, col.names = TRUE, row.names = FALSE)

  svg(paste0(args$prefix, '.variance_explained.svg'), width = 7, height = 7)
  plot(plot_svd_scree(analysis_results$svd_results))
  dev.off()

  svg(paste0(args$prefix, '.svd_significance.svg'), width = 7, height = 7)
  plot(plot_significance(analysis_results$svd_results))
  dev.off()

  svg(paste0(args$prefix, '.sva_significance.svg'), width = 7, height = 7)
  plot(plot_significance(analysis_results$sva_results))
  dev.off()
}
