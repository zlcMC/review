#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: R (epn2_r)
#     language: R
#     name: ir_epn2_r
# ---

# %% [markdown]
# fig03_subclone_undiff_diff_score.R — 对应论文 Figure 3
# 目的：
#   - 在 GTE009 恶性细胞中，基于作者提供的 CNV_Cluster (Subclone_1 / Subclone_2)
#     和 cell_type（EpC-like / NSC-like 等）推出 "分化 / 未分化" 基因签名
#   - 计算 trajectory_score = undiff_score - diff_score
#   - 对全部 4 个样本可视化
# 论文方法摘要 (Methods - Trajectory score)：
#   1) 在 GTE009 恶性细胞里，比较 Subclone_1 vs Subclone_2 的 DEG（FindMarkers / MAST）
#   2) 交集 EpC-like 特异标志 → diff_genes；交集 NSC-like 特异标志 → undiff_genes
#   3) 用 AddModuleScore 给每个细胞计算两套分数
#   4) trajectory_score = undiff - diff
# 运行：conda run -n epn2_r Rscript projectmd/fig03_subclone_undiff_diff_score.R

# %%
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', script_args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
    script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
    setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')

# %%
suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
})
options(Seurat.object.assay.version = 'v3')

# %%
SAMPLES <- c('GTE001','GTE002','GTE009','GTE012')

# %%
# ---- 1) 从 GTE009 推出签名 ----
gte009_path <- output_path('fig1/GTE009_seurat.rds')
stopifnot(file.exists(gte009_path))
obj <- readRDS(gte009_path)

# %%
meta_cols <- colnames(obj@meta.data)
stopifnot('CNV_Cluster' %in% meta_cols, 'Final_cluster' %in% meta_cols)
# 论文 'EpC-like' 对应数据的 'Epe-like'（ependymocyte-like）

# %%
# 仅恶性细胞：作者在 CNV_level 里标出了拷贝数异常的；这里用 Brief_cluster / cell_type 区分
# 论文 Fig3 签名只在恶性细胞 (Subclone_1/2) 间比较
mal <- subset(obj, subset = CNV_Cluster %in% c('Subclone_1','Subclone_2'))
cat(sprintf('GTE009 恶性细胞数: %d\n', ncol(mal)))

# %%
Idents(mal) <- 'CNV_Cluster'
deg <- FindMarkers(mal, ident.1 = 'Subclone_1', ident.2 = 'Subclone_2',
                   only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,
                   test.use = 'wilcox')
deg$gene <- rownames(deg)
write.csv(deg, output_path('fig3/GTE009_subclone_DEG.csv'), row.names = FALSE)

# %%
# Subclone_1 ↑ 对应更分化 (EpC-like)；Subclone_2 ↑ 对应未分化 (NSC-like)
sub1_up <- deg %>% filter(avg_log2FC >  0.5, p_val_adj < 0.05) %>% pull(gene)
sub2_up <- deg %>% filter(avg_log2FC < -0.5, p_val_adj < 0.05) %>% pull(gene)

# %%
# 再用 Final_cluster 交集（Epe-like=EpC-like / NSC-like）
Idents(mal) <- 'Final_cluster'
ct <- unique(as.character(mal$Final_cluster))
cat('Final_cluster 值:', paste(ct, collapse=', '), '\n')

# %%
diff_genes   <- sub1_up
undiff_genes <- sub2_up
epc_label <- grep('^Epe', ct, value=TRUE)
nsc_label <- grep('^NSC', ct, value=TRUE)
if (length(epc_label) && length(nsc_label)) {
    epc_mk <- FindMarkers(mal, ident.1 = epc_label[1],
                          only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    nsc_mk <- FindMarkers(mal, ident.1 = nsc_label[1],
                          only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    diff_genes   <- intersect(sub1_up, rownames(epc_mk)[epc_mk$p_val_adj < 0.05])
    undiff_genes <- intersect(sub2_up, rownames(nsc_mk)[nsc_mk$p_val_adj < 0.05])
}
cat(sprintf('diff_genes (分化/EpC-like): %d 个\n', length(diff_genes)))
cat(sprintf('undiff_genes (未分化/NSC-like): %d 个\n', length(undiff_genes)))

# %%
dir.create(output_path('fig3'), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(diff = diff_genes, undiff = undiff_genes),
        output_path('fig3/trajectory_signatures.rds'))
writeLines(diff_genes,   output_path('fig3/diff_genes.txt'))
writeLines(undiff_genes, output_path('fig3/undiff_genes.txt'))

# %%
# ---- 2) 对每个样本打分并可视化 ----
score_one <- function(s) {
    cat('\n==', s, '==\n')
    rds <- output_path(paste0('fig1/', s, '_seurat.rds'))
    if (!file.exists(rds)) { cat('  缺少', rds, '\n'); return(invisible()) }
    ob <- readRDS(rds)
    ob <- AddModuleScore(ob, features = list(diff_genes),   name = 'diff_score')
    ob <- AddModuleScore(ob, features = list(undiff_genes), name = 'undiff_score')
    ob$trajectory_score <- ob$undiff_score1 - ob$diff_score1

    p <- FeaturePlot(ob, reduction='tsne', features='trajectory_score',
                     cols=c('navy','white','firebrick')) +
         ggtitle(paste(s, 'Trajectory score (undiff - diff)'))
    ggsave(output_path(paste0('fig3/', s, '_trajectory_score.pdf')),
           p, width = 6, height = 5)
    saveRDS(ob@meta.data[,c('diff_score1','undiff_score1','trajectory_score')],
            output_path(paste0('fig3/', s, '_scores.rds')))
    cat('  done\n')
}

# %%
for (s in SAMPLES) score_one(s)
cat('\n✓ Figure 3 trajectory score 完成\n')
