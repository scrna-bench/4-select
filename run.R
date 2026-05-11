#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(HDF5Array)
})

# Parse command line arguments
parser <- ArgumentParser(description="OmniBenchmark module")

# Required by OmniBenchmark
parser$add_argument("--output_dir", dest="output_dir", type="character", required=TRUE,
                   help="Output directory for results")
parser$add_argument("--name", dest="name", type="character", required=TRUE,
                   help="Module name/identifier")

# Stage-specific inputs
parser$add_argument("--normalized.h5", dest="input_h5",
                   type="character", nargs="+", required=TRUE,
                   help="Input: normalized.h5")
parser$add_argument("--rawdata.h5ad", dest="rawdata_h5ad",
                   type="character", required=FALSE,
                   help="Raw h5ad file (used for raw counts by seurat_vst methods)")
parser$add_argument("--filtered.cellids", dest="filtered_cellids",
                   type="character", required=FALSE,
                   help="Filtered cell IDs (gzipped text)")

parser$add_argument("--selection_type", dest="selection_type",
                    type="character", help="Gene selection method")
parser$add_argument("--number_selected", dest="number_selected",
                    type="character", help="Number of genes to select")
parser$add_argument("--batch_variable", dest="batch_variable",
                   type="character", default=NULL,
                   help="colData column name to use as batch variable")

args <- parser$parse_args()

cat("Full command: ", paste0(commandArgs(), collapse = " "), "\n")
cat("output_dir:", args$output_dir, "\n")
cat("name:", args$name, "\n")
cat("selection_type:", args$selection_type, "\n")
args$number_selected <- as.integer(args$number_selected)
cat("number_selected:", args$number_selected, "\n")
cat("input_h5:", args$input_h5, "\n")
cat("rawdata_h5ad:", args$rawdata_h5ad, "\n")
cat("filtered_cellids:", args$filtered_cellids, "\n")
cat("batch_variable:", args$batch_variable, "\n")

if (args$selection_type == "seurat_vst") {
  require(Seurat)
  so <- read_h5ad(args$rawdata_h5ad, as = "Seurat")
  cellids <- readLines(gzfile(args$filtered_cellids))
  so <- subset(so, cells = cellids)
  cat("dim(so) after filtering:", dim(so), "\n")
  so <- FindVariableFeatures(so, selection.method = "vst",
        nfeatures = args$number_selected)
  sel_feats <- VariableFeatures(so)
} else if (args$selection_type == "seurat_vst_batch") {
  require(Seurat)
  so <- read_h5ad(args$rawdata_h5ad, as = "Seurat")
  cellids <- readLines(gzfile(args$filtered_cellids))
  so <- subset(so, cells = cellids)
  cat("dim(so) after filtering:", dim(so), "\n")
  batch_col <- args$batch_variable
  batches <- unique(so[[batch_col, drop = TRUE]])
  cat("batches:", batches, "\n")
  seurat_list <- lapply(batches, function(b) {
    cells_b <- colnames(so)[so[[batch_col, drop = TRUE]] == b]
    sub_so <- subset(so, cells = cells_b)
    FindVariableFeatures(sub_so, selection.method = "vst",
                         nfeatures = args$number_selected)
  })
  sel_feats <- SelectIntegrationFeatures(seurat_list,
                                         nfeatures = args$number_selected)
} else {
  stop("incorrect 'selection_type' specified")
}

cat("length(sel_feats):", length(sel_feats), "\n")

m <- TENxMatrix(args$input_h5, group = "matrix")
m <- as(m, "dgCMatrix")

output_file <- file.path(args$output_dir, paste0(args$name, "_normalized_selected.h5"))
cat("output_file:", output_file, "\n")
writeTENxMatrix(m[sel_feats,], output_file, group="matrix")
file.info(output_file)[,c("size", "ctime")]

