# Main functions for the OmniBenchmark module

#' Process data using parsed command-line arguments
#'
#' @param args Parsed arguments containing:
#'   - output_dir: Output directory path
#'   - name: Module name
#'   - normalized_h5: Input files for normalized.h5 (CLI: --normalized.h5)
#'   - filtered_cellids: Input files for filtered.cellids (CLI: --filtered.cellids)
#'
#' @note Input IDs with dots (e.g., 'data.raw') are converted to underscores
#'   in R variable names (e.g., 'data_raw') but preserve dots in CLI args.
process_data <- function(args) {
  # Create output directory if needed
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Processing module:", args$name, "\n")

  # Access stage inputs
  normalized_h5_files <- args$normalized_h5
  cat("  normalized.h5:", normalized_h5_files, "\n")
  filtered_cellids_files <- args$filtered_cellids
  cat("  filtered.cellids:", filtered_cellids_files, "\n")

  # TODO: Implement your processing logic here
  # Example: Read inputs, process, write outputs

  # Write a simple output file
  output_file <- file.path(args$output_dir, paste0(args$name, "_result.txt"))
  output_lines <- c(
    paste("Processed module:", args$name),
    paste("normalized.h5:", length(normalized_h5_files), "file(s)"),
    paste("filtered.cellids:", length(filtered_cellids_files), "file(s)")
  )
  writeLines(output_lines, output_file)

  cat("Results written to:", output_file, "\n")
}
