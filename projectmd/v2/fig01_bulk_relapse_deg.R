#!/usr/bin/env Rscript
# fig01_bulk_relapse_deg.R (v2) — Figure 1F: bulk EPN relapse vs primary DEG vs sc DEG.
# Uses GSE64415 (Pajtler 2015, 209 EPN, primary/relapse labels).
#
# v2 changes:
#   * Uses helpers::init_workspace + helpers::load_gse_expression to drop the
#     duplicated GEO loading boilerplate (~30 lines).
#   * Writes outputs to output/v2/fig1_bulk/ instead of output/fig1_bulk/.
#   * Reads sc DEG from output/v2/fig1_supp/ first, falling back to output/fig1_supp/.

this_script <- sub('^--file=', '',
                   grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(this_script), 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(VennDiagram); library(grid); library(fgsea); library(limma)
})

out <- v2_dir('fig1_bulk')

# ---- 1) load GSE64415 ----
dat <- load_gse_expression(
  matrix_path = raw_data_path('bulk/GSE64415_series_matrix.txt.gz'),
  extra_columns = list(
    primary_relapse = 'primary/relapse',
    subgroup        = 'molecular subgroup'
  ),
  sample_filter_key    = 'primary_relapse',
  sample_filter_values = c('primary', 'relapse')
)
expr <- dat$expr
group <- factor(dat$meta$primary_relapse, levels = c('primary', 'relapse'))
cat(sprintf('expr: %d × %d (primary=%d, relapse=%d)\n',
            nrow(expr), ncol(expr), sum(group == 'primary'), sum(group == 'relapse')))

# ---- 2) limma DEG ----
design <- model.matrix(~ group)
fit    <- eBayes(lmFit(expr, design))
deg    <- topTable(fit, coef = 2, number = Inf, adjust.method = 'BH')
deg$gene <- rownames(deg)
write.csv(deg, file.path(out, 'bulk_DEG_relapse_vs_primary.csv'), row.names = FALSE)
cat('Top 10 bulk DEG (relapse vs primary):\n')
print(head(deg[, c('gene', 'logFC', 'P.Value', 'adj.P.Val')], 10))

# ---- 3) sc DEG intersection ----
sc_path <- dep_file('fig1_supp', 'GTE009_Sub1vs2_DEG.csv')
if (!file.exists(sc_path)) {
  cands <- list.files(dirname(sc_path), pattern = 'DEG|marker', full.names = TRUE)
  if (length(cands)) sc_path <- cands[1]
}
if (file.exists(sc_path)) {
  sc <- read.csv(sc_path)
  cat(sprintf('sc DEG: %s (%d rows)\n', sc_path, nrow(sc)))
  fc_col <- grep('log2FC|logFC|avg_log', colnames(sc), value = TRUE, ignore.case = TRUE)[1]
  p_col  <- grep('p_val_adj|adj.P|padj',  colnames(sc), value = TRUE, ignore.case = TRUE)[1]
  g_col  <- intersect(c('gene', 'symbol', 'X', 'Gene'), colnames(sc))[1]
  if (is.na(g_col)) g_col <- colnames(sc)[1]

  sc_up   <- sc[[g_col]][sc[[fc_col]] >  0.25 & sc[[p_col]] < 0.05]
  sc_dn   <- sc[[g_col]][sc[[fc_col]] < -0.25 & sc[[p_col]] < 0.05]
  bulk_up <- deg$gene[deg$logFC >  0.5 & deg$adj.P.Val < 0.1]
  bulk_dn <- deg$gene[deg$logFC < -0.5 & deg$adj.P.Val < 0.1]
  cat(sprintf('sc up=%d dn=%d  | bulk up=%d dn=%d\n',
              length(sc_up), length(sc_dn), length(bulk_up), length(bulk_dn)))

  ov_up <- intersect(sc_up, bulk_up); ov_dn <- intersect(sc_dn, bulk_dn)
  write.csv(data.frame(gene = ov_up), file.path(out, 'sc_vs_bulk_overlap_up.csv'), row.names = FALSE)
  write.csv(data.frame(gene = ov_dn), file.path(out, 'sc_vs_bulk_overlap_dn.csv'), row.names = FALSE)
  cat(sprintf('Overlap up=%d dn=%d\n', length(ov_up), length(ov_dn)))

  for (sfx in c('up', 'dn')) {
    pdf(file.path(out, sprintf('Fig1F_venn_%s.pdf', sfx)), width = 5, height = 5)
    fill <- if (sfx == 'up') c('#E64B35', '#4DBBD5') else c('#3C5488', '#00A087')
    sets <- if (sfx == 'up') list(sc_Sub1_up = sc_up, bulk_relapse_up = bulk_up)
            else            list(sc_Sub1_dn = sc_dn, bulk_relapse_dn = bulk_dn)
    grid.draw(venn.diagram(sets, filename = NULL, fill = fill, alpha = 0.5))
    dev.off()
  }

  ranks <- deg$logFC; names(ranks) <- deg$gene
  ranks <- sort(ranks[!is.na(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
  pw <- list(sc_Sub1_up = sc_up, sc_Sub1_dn = sc_dn)
  pw <- pw[lengths(pw) >= 5]
  if (length(pw)) {
    res <- fgsea(pathways = pw, stats = ranks, minSize = 5)
    res_out <- as.data.frame(res)
    res_out$leadingEdge <- sapply(res_out$leadingEdge, paste, collapse = ';')
    write.csv(res_out, file.path(out, 'GSEA_sc_in_bulk.csv'), row.names = FALSE)
    print(res[, c('pathway', 'pval', 'padj', 'NES', 'size')])
  }
} else {
  cat('No sc DEG found, skipping intersection analysis.\n')
}

cat('\n[v2] Fig 1F bulk DEG done →', out, '\n')
