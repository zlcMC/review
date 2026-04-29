#!/usr/bin/env Rscript
# v2 of fig05_cellchat_downstream.R — entry that runs the two figure scripts.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

message('Step 1/2: full-information CellChat figures -> ', v2_dir('fig5_fullsize'))
source(file.path(this_dir, 'fig05_cellchat_fullsize_figures.R'), local = FALSE)

message('Step 2/2: reader-friendly CellChat figures -> ', v2_dir('fig5_readable'))
source(file.path(this_dir, 'fig05_cellchat_readable_figures.R'), local = FALSE)

message('Done.')
