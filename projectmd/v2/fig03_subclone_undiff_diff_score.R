#!/usr/bin/env Rscript
# v2 of fig03_subclone_undiff_diff_score.R
# Build undiff/diff trajectory signatures from GTE009 and score all 4 samples.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr)
})
options(Seurat.object.assay.version = 'v3')

SAMPLES <- c('GTE001', 'GTE002', 'GTE009', 'GTE012')

# ---- 1) Derive signatures from GTE009 ----
obj <- readRDS(dep_file('fig1', 'GTE009_seurat.rds'))
stopifnot(all(c('CNV_Cluster', 'Final_cluster') %in% colnames(obj@meta.data)))

mal <- subset(obj, subset = CNV_Cluster %in% c('Subclone_1', 'Subclone_2'))
cat(sprintf('GTE009 malignant cells: %d\n', ncol(mal)))

Idents(mal) <- 'CNV_Cluster'
deg <- FindMarkers(mal, ident.1 = 'Subclone_1', ident.2 = 'Subclone_2',
                   only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
                   test.use = 'wilcox')
deg$gene <- rownames(deg)
write.csv(deg, v2_file('fig3', 'GTE009_subclone_DEG.csv'), row.names = FALSE)

sub1_up <- deg %>% filter(avg_log2FC >  0.5, p_val_adj < 0.05) %>% pull(gene)
sub2_up <- deg %>% filter(avg_log2FC < -0.5, p_val_adj < 0.05) %>% pull(gene)

Idents(mal) <- 'Final_cluster'
ct <- unique(as.character(mal$Final_cluster))
cat('Final_cluster values:', paste(ct, collapse = ', '), '\n')

epc <- grep('^Epe', ct, value = TRUE)
nsc <- grep('^NSC', ct, value = TRUE)
diff_genes   <- sub1_up
undiff_genes <- sub2_up
if (length(epc) && length(nsc)) {
  epc_mk <- FindMarkers(mal, ident.1 = epc[1], only.pos = TRUE,
                        min.pct = 0.25, logfc.threshold = 0.25)
  nsc_mk <- FindMarkers(mal, ident.1 = nsc[1], only.pos = TRUE,
                        min.pct = 0.25, logfc.threshold = 0.25)
  diff_genes   <- intersect(sub1_up, rownames(epc_mk)[epc_mk$p_val_adj < 0.05])
  undiff_genes <- intersect(sub2_up, rownames(nsc_mk)[nsc_mk$p_val_adj < 0.05])
}
cat(sprintf('diff (EpC-like): %d   undiff (NSC-like): %d\n',
            length(diff_genes), length(undiff_genes)))

saveRDS(list(diff = diff_genes, undiff = undiff_genes),
        v2_file('fig3', 'trajectory_signatures.rds'))
writeLines(diff_genes,   v2_file('fig3', 'diff_genes.txt'))
writeLines(undiff_genes, v2_file('fig3', 'undiff_genes.txt'))

# ---- 2) Score every sample ----
score_one <- function(s) {
  cat('\n==', s, '==\n')
  rds <- dep_file('fig1', paste0(s, '_seurat.rds'))
  if (!file.exists(rds)) { cat('  missing', rds, '\n'); return(invisible()) }
  ob <- readRDS(rds)
  ob <- AddModuleScore(ob, features = list(diff_genes),   name = 'diff_score')
  ob <- AddModuleScore(ob, features = list(undiff_genes), name = 'undiff_score')
  ob$trajectory_score <- ob$undiff_score1 - ob$diff_score1

  ggsave(v2_file('fig3', paste0(s, '_trajectory_score.pdf')),
         FeaturePlot(ob, reduction = 'tsne', features = 'trajectory_score',
                     cols = c('navy', 'white', 'firebrick')) +
           ggtitle(paste(s, 'Trajectory score (undiff - diff)')),
         width = 6, height = 5)
  saveRDS(ob@meta.data[, c('diff_score1', 'undiff_score1', 'trajectory_score')],
          v2_file('fig3', paste0(s, '_scores.rds')))
  cat('  done\n')
}

for (s in SAMPLES) score_one(s)
cat('\n✓ Figure 3 trajectory score (v2) done\n')
