#!/usr/bin/env Rscript
# fig05_cellchat_downstream.R
# Clean downstream entry point for Fig 5 / Supp Fig 8 CellChat outputs.
#
# Input:
#   output/fig5/cellchat_merged.rds from fig05_cellchat_prepare.R
# Outputs:
#   output/fig5_fullsize/  full-information mainline figures
#   output/fig5_readable/  reader-friendly interpretation figures
#
# Run:
#   conda run -n epn2_r Rscript projectmd/fig05_cellchat_downstream.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
    script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
    setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')

message('Step 1/2: writing full-information CellChat figures to output/fig5_fullsize')
source(code_path('fig05_cellchat_fullsize_figures.R'), local = FALSE)

message('Step 2/2: writing reader-friendly CellChat figures to output/fig5_readable')
source(code_path('fig05_cellchat_readable_figures.R'), local = FALSE)

message('Done. Final CellChat downstream outputs:')
message('  - ', output_path('fig5_fullsize'))
message('  - ', output_path('fig5_readable'))
