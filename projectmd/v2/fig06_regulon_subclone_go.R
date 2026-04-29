#!/usr/bin/env Rscript
# task 2 — Subclone-上调 regulon 的 TF symbol 跑 GO/BP 富化，给 Fig 6 增量加生物学解释。

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

out <- v2_dir('fig6')

clean_tf <- function(x) sub('\\(\\+\\)$', '', x)

run_go <- function(csv_in, label, fdr_thr = 0.05, delta_col, sign = 'up') {
  df <- read.csv(csv_in, check.names = FALSE)
  if (sign == 'up') sub <- df[df$fdr < fdr_thr & df[[delta_col]] > 0, ]
  else              sub <- df[df$fdr < fdr_thr & df[[delta_col]] < 0, ]
  tfs <- unique(clean_tf(sub$regulon))
  cat(sprintf('[%s] %d significant regulons (FDR<%.2f, sign=%s)\n',
              label, length(tfs), fdr_thr, sign))
  if (length(tfs) < 5) { cat('Too few TFs, skip GO\n'); return(invisible()) }
  enrich_go_dotplot(
    genes      = tfs,
    title      = sprintf('GO BP — %s', label),
    out_pdf    = file.path(out, sprintf('Fig6_regulon_%s_GO_BP.pdf', label)),
    out_png    = file.path(out, sprintf('Fig6_regulon_%s_GO_BP.png', label)),
    csv_path   = file.path(out, sprintf('Fig6_regulon_%s_GO_BP.csv', label)),
    genes_txt  = file.path(out, sprintf('Fig6_regulon_%s_TFs.txt', label)),
    pcutoff    = 0.05, qcutoff = 0.1, csv_writer = 'base'
  )
}

run_go(file.path(out, 'regulon_subclone_diff_GTE009.csv'),
       'GTE009_Subclone1_up', delta_col = 'delta_Subclone_1_minus_Subclone_2', sign = 'up')
run_go(file.path(out, 'regulon_subclone_diff_GTE009.csv'),
       'GTE009_Subclone2_up', delta_col = 'delta_Subclone_1_minus_Subclone_2', sign = 'down')
run_go(file.path(out, 'regulon_subclone_diff_pooled.csv'),
       'Pooled_Subclone1_up', delta_col = 'delta_Subclone_1_minus_Subclone_2', sign = 'up')

cat('✓ Subclone × regulon GO done\n')
