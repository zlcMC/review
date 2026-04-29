#!/usr/bin/env Rscript
# fig02_slingshot_basic.R (v2) — Sanity-check Slingshot on GTE009 malignant cells
# 输出 pseudotime + 4-panel 概览图，供后续 fig02c_slingshot_final_panel.R 复核。
# 运行：conda run -n epn2_r Rscript projectmd/v2/fig02_slingshot_basic.R

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(SingleCellExperiment)
  library(ggplot2); library(patchwork)
})

v2_dir('fig2')
obj <- readRDS(dep_file('fig1', 'GTE009_seurat.rds'))
mal <- subset(obj, subset = CNV_Cluster %in% c('Subclone_1', 'Subclone_2'))
cat('GTE009 malignant cells:', ncol(mal), '\n')

sce <- as.SingleCellExperiment(mal, assay = 'RNA')
sds <- slingshot(sce, reducedDim = 'PCA',
                 clusterLabels = mal$Final_cluster, start.clus = 'NSC-like')
saveRDS(sds, v2_file('fig2', 'GTE009_slingshot.rds'))

mal$pseudotime <- slingPseudotime(sds)[, 1]

p1 <- FeaturePlot(mal, reduction='tsne', features='pseudotime',
                  cols=c('navy','white','firebrick')) +
      ggtitle('Fig 2C - Slingshot pseudotime (NSC-like start)')
p2 <- DimPlot(mal, reduction='tsne', group.by='Final_cluster', label=TRUE) + ggtitle('Final_cluster')
p3 <- DimPlot(mal, reduction='tsne', group.by='CNV_Cluster',  label=TRUE) + ggtitle('CNV_Cluster')
p4 <- if ('CytoTRACE' %in% colnames(mal@meta.data)) {
  FeaturePlot(mal, reduction='tsne', features='CytoTRACE',
              cols=c('blue','white','red')) + ggtitle('CytoTRACE')
} else patchwork::plot_spacer()

ggsave(v2_file('fig2', 'Fig2C_GTE009_slingshot.pdf'), (p1 | p2) / (p3 | p4),
       width = 14, height = 10)
saveRDS(mal@meta.data[, c('Final_cluster','CNV_Cluster','pseudotime','CytoTRACE')],
        v2_file('fig2', 'GTE009_pseudotime.rds'))
cat('✓ Fig 2C basic Slingshot 完成\n')
