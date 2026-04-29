#!/usr/bin/env Rscript
# v2 of fig04_gojo_scrna_deg_go.R
# GO enrichment for Gojo GSE141460 PF-A recurrent-high malignant pseudobulk DEGs.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({ library(readr); library(dplyr) })

deg_path <- dep_file('fig4_gojo_scrna',
                     'GSE141460_malignant_pseudobulk_DEG_PFA_only.csv')
if (!file.exists(deg_path)) {
  stop('Missing DEG file. Run fig04_gojo_scrna_validation.py first: ', deg_path)
}

deg <- read_csv(deg_path, show_col_types = FALSE) %>%
  filter(!is.na(ttest_p), recurrent_minus_primary > 0)

candidate <- deg %>% filter(ttest_p < 0.05, recurrent_minus_primary > 0.25)
if (nrow(candidate) < 10) {
  candidate <- deg %>% arrange(ttest_p, desc(recurrent_minus_primary)) %>% head(150)
}

enrich_go_dotplot(
  genes     = candidate$gene,
  title     = 'Gojo GSE141460 PF-A recurrent-high pseudobulk GO',
  out_pdf   = v2_file('fig4_gojo_scrna', 'Fig4D_GSE141460_PFA_recurrent_high_GO_BP_scrna.pdf'),
  out_png   = v2_file('fig4_gojo_scrna', 'Fig4D_GSE141460_PFA_recurrent_high_GO_BP_scrna.png'),
  csv_path  = v2_file('fig4_gojo_scrna', 'GSE141460_PFA_recurrent_high_GO_BP.csv'),
  genes_txt = v2_file('fig4_gojo_scrna', 'GSE141460_recurrent_high_GO_input_genes.txt')
)
