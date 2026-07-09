#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(HDF5Array)
  library(anndataR)
  library(SingleCellExperiment)
  library(scry)
})

run_select <- function(args) {
  cellids <- readLines(gzfile(args$filtered_cellids))

  if (args$selection_type == "scry_deviance") {
    sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment")
    sce <- sce[, cellids]
    cat(sprintf("  dim(sce) after filtering: %d x %d\n", nrow(sce), ncol(sce)))
    sce <- devianceFeatureSelection(sce, assay = "counts", nkeep = args$number_selected)
    rownames(sce)
  } else {
    stop("Unsupported selection_type: ", args$selection_type)
  }
}

main <- function() {
  parser <- ArgumentParser(description = "4-select: scry-based feature selection")
  parser$add_argument("--output_dir",       type = "character", required = TRUE, help = "output directory")
  parser$add_argument("--name",             type = "character", required = TRUE, help = "module name")
  parser$add_argument("--normalized_h5",    type = "character", dest = "input_h5",     nargs = "+", help = "normalized HDF5 matrix")
  parser$add_argument("--rawdata_h5ad",     type = "character", nargs = "+",           help = "raw data h5ad")
  parser$add_argument("--filtered_cellids", type = "character", help = "filtered cell IDs gz")
  parser$add_argument("--properties_info",  type = "character", help = "properties yaml (unused)")
  parser$add_argument("--selection_type",   type = "character", help = "selection method")
  parser$add_argument("--number_selected",  type = "integer",   help = "number of features to select")
  args <- parser$parse_args()

  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "input_h5", "rawdata_h5ad",
              "filtered_cellids", "selection_type", "number_selected")) {
    cat(sprintf("  %s: %s\n", k, args[[k]]))
  }

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  sel_feats <- run_select(args)
  cat("length(sel_feats):", length(sel_feats), "\n")

  m <- TENxMatrix(args$input_h5, group = "matrix")
  m <- as(m, "dgCMatrix")

  out <- file.path(args$output_dir, paste0(args$name, "_normalized_selected.h5"))
  cat("output_file:", out, "\n")
  writeTENxMatrix(m[sel_feats, ], out, group = "matrix")
  cat(sprintf("  wrote: %s\n", out))
  print(file.info(out)[, c("size", "ctime")])
}

if (sys.nframe() == 0L) {
  main()
}
