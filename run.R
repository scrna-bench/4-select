#!/usr/bin/env Rscript

library(argparse)
library(HDF5Array)

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

parser$add_argument("--selection_type", dest="selection_type", 
                    type="character", help="Input file")
parser$add_argument("--number_selected", dest="number_selected", 
                    type="character", help="Input file")

args <- parser$parse_args()

cat("Full command: ", paste0(commandArgs(), collapse = " "), "\n")
cat("output_dir:", args$output_dir, "\n")
cat("name:", args$name, "\n")
cat("selection_type:", args$selection_type, "\n")
args$number_selected <- as.integer(args$number_selected)
cat("number_selected:", args$number_selected, "\n")
cat("input_h5:", args$input_h5, "\n")
cat("cellids:", args$cellids, "\n")

cellids <- readLines(gzfile(args$cellids))
cat("length(cellids):", length(cellids), "\n")

m <- TENxMatrix(args$input_h5, group = "matrix")
m <- as(m, "dgCMatrix") # read into memory
cat("dim(m):", dim(m), "\n")

if (args$selection_type == "seurat_vst") {
  require(Seurat)
  so <- CreateSeuratObject(counts = m, assay = "normalized")
  so <- FindVariableFeatures(so, selection.method = "vst", 
        nfeatures = args$number_selected)
  sel_feats <- VariableFeatures(so)
} else if (args$selection_type == "scrapper_modelGeneVariances") {
  #require(SingleCellExperiment)
  require(scrapper)
  #sce <- SingleCellExperiment(list(normalized=m)) # not needed, actually
  gene.var <- modelGeneVariances(m)
  hvgs <- chooseHighlyVariableGenes(gene.var$statistics$residuals,
                                    top = args$number_selected)
  sel_feats <- rownames(m)[hvgs]
} else {
  errorCondition("incorrect 'selection_type' specified")
}

cat("length(sel_feats):", length(sel_feats), "\n")

output_file <- file.path(args$output_dir, paste0(args$name, "_selected.txt.gz"))
writeLines(sel_feats, gzfile(output_file))

