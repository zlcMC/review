#!/usr/bin/env Rscript
# v2 of fig04_gillen_scrna_deg_go.R
# GO enrichment for Gillen GSE125969 recurrent-high neoplastic pseudobulk DEGs.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({ library(readr); library(dplyr) })

deg_path <- dep_file('fig4_scrna',
                     'GSE125969_neoplastic_pseudobulk_DEG_all_subtypes.csv')
if (!file.exists(deg_path)) {
  stop('Missing DEG file. Run fig04_gillen_scrna_pseudobulk_lr.py first: ', deg_path)
}

deg <- read_csv(deg_path, show_col_types = FALSE) %>%
  filter(!is.na(ttest_p), recurrent_minus_primary_logFC > 0)

candidate <- deg %>% filter(ttest_p < 0.05, recurrent_minus_primary_logFC > 0.25)
if (nrow(candidate) < 10) {
  candidate <- deg %>% arrange(ttest_p, desc(recurrent_minus_primary_logFC)) %>% head(150)
}

enrich_go_dotplot(
  genes     = candidate$gene,
  title     = 'GSE125969 recurrent-high neoplastic pseudobulk GO',
  out_pdf   = v2_file('fig4_scrna', 'Fig4D_GSE125969_recurrent_high_GO_BP_scrna.pdf'),
  out_png   = v2_file('fig4_scrna', 'Fig4D_GSE125969_recurrent_high_GO_BP_scrna.png'),
  csv_path  = v2_file('fig4_scrna', 'GSE125969_recurrent_high_GO_BP.csv'),
  genes_txt = v2_file('fig4_scrna', 'GSE125969_recurrent_high_GO_input_genes.txt')
)
