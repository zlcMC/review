#!/usr/bin/env Rscript
# fig02_slingshot_basic.R — 替代 monocle2 (dplyr 不兼容)
# 论文 Fig 2C: GTE009 恶性细胞的差化轨迹
# 改用 slingshot (与 monocle DDRTree 同属 trajectory inference)
#   - 在 GTE009 恶性细胞 (Subclone_1/2) 上跑
#   - 用作者的 PCA + Final_cluster 作为聚类
#   - slingshot 输出 pseudotime + 主曲线
# 运行：conda run -n epn2_r Rscript projectmd/fig02_slingshot_basic.R

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', script_args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
      script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
      setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')

suppressPackageStartupMessages({
    library(Seurat); library(slingshot); library(SingleCellExperiment)
    library(ggplot2); library(dplyr); library(patchwork)
})

dir.create(output_path('fig2'), showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(output_path('fig1/GTE009_seurat.rds'))
mal <- subset(obj, subset = CNV_Cluster %in% c('Subclone_1','Subclone_2'))
cat('GTE009 恶性细胞:', ncol(mal), '\n')

# 转 SingleCellExperiment
sce <- as.SingleCellExperiment(mal, assay = 'RNA')
# slingshot 在 PCA + Final_cluster 上做轨迹
sds <- slingshot(sce, reducedDim = 'PCA',
                 clusterLabels = mal$Final_cluster,
                 start.clus = 'NSC-like')   # 论文：NSC-like 起源
saveRDS(sds, output_path('fig2/GTE009_slingshot.rds'))

# pseudotime 取第一条主曲线
pt <- slingPseudotime(sds)
mal$pseudotime <- pt[, 1]

# tSNE 上画 pseudotime + Final_cluster + Subclone
p1 <- FeaturePlot(mal, reduction='tsne', features='pseudotime',
                  cols=c('navy','white','firebrick')) +
      ggtitle('Fig 2C - Slingshot pseudotime (NSC-like start)')
p2 <- DimPlot(mal, reduction='tsne', group.by='Final_cluster', label=TRUE) +
      ggtitle('Final_cluster')
p3 <- DimPlot(mal, reduction='tsne', group.by='CNV_Cluster', label=TRUE) +
      ggtitle('CNV_Cluster (Subclone)')
p4 <- if ('CytoTRACE' %in% colnames(mal@meta.data)) {
    FeaturePlot(mal, reduction='tsne', features='CytoTRACE',
                cols=c('blue','white','red')) + ggtitle('CytoTRACE')
} else { patchwork::plot_spacer() }
g <- (p1 | p2) / (p3 | p4)
ggsave(output_path('fig2/Fig2C_GTE009_slingshot.pdf'), g, width = 14, height = 10)

# 保存 pseudotime 给后续 Figure 2/4 复核步骤
saveRDS(mal@meta.data[, c('Final_cluster','CNV_Cluster','pseudotime','CytoTRACE')],
        output_path('fig2/GTE009_pseudotime.rds'))
cat('✓ Fig 2C (slingshot 替代 monocle) 完成\n')
