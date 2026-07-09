#!/usr/bin/env Rscript
# scry deviance-based feature selection module for omnibenchmark.
#
# Reads raw counts from rawdata.h5ad, computes binomial deviance per gene
# via scry::devianceFeatureSelection, and outputs the top N genes as a
# subsetted normalized HDF5 matrix.

suppressPackageStartupMessages({
  library(argparser)
  library(HDF5Array)
  library(anndataR)
  library(SingleCellExperiment)
  library(scry)
})

# arg parsing
p <- arg_parser("FEAT module — scry deviance-based feature selection")
p <- add_argument(p, "--output_dir",       help = "output directory")
p <- add_argument(p, "--name",             help = "module name")
p <- add_argument(p, "--rawdata_h5ad",     help = "raw h5ad (counts in X)", nargs = "+")
p <- add_argument(p, "--normalized_h5",    help = "normalized HDF5 matrix", nargs = "+")
p <- add_argument(p, "--filtered_cellids", help = "filtered cell IDs gz")
p <- add_argument(p, "--properties_info",  help = "properties yaml (unused)")
p <- add_argument(p, "--selection_type",   help = "selection method")
p <- add_argument(p, "--number_selected",  type = "integer", help = "number of features to select")
args <- parse_args(p)

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in seq_along(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  cellids <- readLines(gzfile(args$filtered_cellids))

  if (args$selection_type == "scry_deviance") {
    sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment")
    sce <- sce[, cellids]
    cat(sprintf("  dim(sce) after filtering: %d x %d\n", nrow(sce), ncol(sce)))
    sce <- devianceFeatureSelection(sce, assay = "counts", nkeep = args$number_selected)
    sel_feats <- rownames(sce)
  } else {
    stop("Unsupported selection_type: ", args$selection_type)
  }

  cat(sprintf("  length(sel_feats): %d\n", length(sel_feats)))

  m <- TENxMatrix(args$normalized_h5, group = "matrix")
  m <- as(m, "dgCMatrix")

  out <- file.path(args$output_dir, paste0(args$name, "_normalized_selected.h5"))
  cat(sprintf("output_file: %s\n", out))
  writeTENxMatrix(m[sel_feats, ], out, group = "matrix")
  cat(sprintf("  wrote: %s\n", out))
  print(file.info(out)[, c("size", "ctime")])
}

if (sys.nframe() == 0L) {
  main()
}
