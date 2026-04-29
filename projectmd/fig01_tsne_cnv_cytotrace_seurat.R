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
# fig01_tsne_cnv_cytotrace_seurat.R — Figure 1 tSNE/CNV/CytoTRACE/Seurat export
# 目的：
#   - 按论文 Methods 读入 GTE001/002/009/012 的 counts + 作者 metadata
#   - 归一化 + HVG(5000) + CellCycle 回归 + PCA(50) + tSNE（论文用 tSNE 不是 UMAP）
#   - 应用作者注释 (Brief_cluster / CNV_Cluster / CNV_level / CytoTRACE)
#   - 输出 Figure 1A-E 对应的 tSNE 图到 output/fig1/
# 运行：conda run -n epn2_r Rscript projectmd/fig01_tsne_cnv_cytotrace_seurat.R

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
    library(patchwork)
})
options(Seurat.object.assay.version = 'v3')

# %%
SAMPLES <- c('GTE001', 'GTE002', 'GTE009', 'GTE012')
GSM <- c(GTE001='GSM5710226', GTE002='GSM5710227',
         GTE009='GSM5710228', GTE012='GSM5710229')

# %%
process_one <- function(s) {
    cat('\n====', s, '====\n')
    cache <- output_path(paste0('fig1/', s, '_seurat.rds'))
    if (file.exists(cache)) {
        cat('  (cached) 载入', cache, '\n')
        return(readRDS(cache))
    }
    t0 <- Sys.time()
    counts_file <- raw_data_path(paste0(GSM[s], '_', s, '_counts.csv.gz'))
    meta_file   <- raw_data_path(paste0(GSM[s], '_', s, '_metadata.csv.gz'))
    df   <- read.csv(counts_file, row.names = 1, check.names = FALSE)
    meta <- read.csv(meta_file,   row.names = 1, check.names = FALSE)
    cat(sprintf('  读入: %d × %d  (%.1fs)\n', nrow(df), ncol(df),
                as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

    obj <- CreateSeuratObject(counts = df, meta.data = meta, min.cells = 10)
    rm(df); invisible(gc())

    # 论文阈值：nFeature ≥ 1500，pct_mt < 12%
    # 注意：作者上传的细胞是已经做过 QC 的，这里复核一遍；metadata 列名是 percent.mt
    if (!'percent.mt' %in% colnames(obj@meta.data)) {
        obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = '^MT-')
    }
    obj <- subset(obj, subset = nFeature_RNA >= 1500 & percent.mt < 12)

    # DoubletFinder 结果已在 meta 的 DF_hi.lo 中
    if ('DF_hi.lo' %in% colnames(obj@meta.data)) {
        obj <- subset(obj, subset = DF_hi.lo == 'Singlet')
    }
    cat(sprintf('  QC+Singlet 后: %d 细胞\n', ncol(obj)))

    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 5000, verbose = FALSE)
    # S.Score / G2M.Score 已在 meta，不再重算；若缺失则补
    if (!all(c('S.Score','G2M.Score') %in% colnames(obj@meta.data))) {
        obj <- CellCycleScoring(obj,
            s.features = cc.genes.updated.2019$s.genes,
            g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
    }
    obj <- ScaleData(obj, vars.to.regress = c('S.Score','G2M.Score'), verbose = FALSE)
    obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
    obj <- RunTSNE(obj, dims = 1:50, seed.use = 42)
    obj <- FindNeighbors(obj, dims = 1:50, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
    cat(sprintf('  PCA+tSNE+cluster 完成  (总耗时 %.1fs)\n',
                as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

    saveRDS(obj, cache)
    cat('  已缓存到', cache, '\n')
    obj
}

# %%
plots_for_sample <- function(obj, s) {
    dir.create(output_path('fig1'), showWarnings = FALSE, recursive = TRUE)
    # Figure 1 面板（tSNE 投影）
    p_bc  <- DimPlot(obj, reduction='tsne', group.by='Brief_cluster',
                     label=TRUE, pt.size=0.3) + ggtitle(paste(s, 'Brief_cluster'))
    p_cnv <- if ('CNV_level' %in% colnames(obj@meta.data))
        FeaturePlot(obj, reduction='tsne', features='CNV_level',
                    cols=c('grey90','firebrick')) + ggtitle(paste(s,'CNV level (Fig 1A)'))
        else patchwork::plot_spacer()
    p_cyt <- if ('CytoTRACE' %in% colnames(obj@meta.data))
        FeaturePlot(obj, reduction='tsne', features='CytoTRACE',
                    cols=c('blue','white','red')) + ggtitle(paste(s,'CytoTRACE (Fig 1B)'))
        else patchwork::plot_spacer()
    p_sub <- if ('CNV_Cluster' %in% colnames(obj@meta.data))
        DimPlot(obj, reduction='tsne', group.by='CNV_Cluster',
                label=TRUE, pt.size=0.3) + ggtitle(paste(s,'Subclone (Fig 1C/D)'))
        else patchwork::plot_spacer()
    g <- (p_bc | p_cnv) / (p_cyt | p_sub)
    ggsave(output_path(paste0('fig1/', s, '_tsne_panels.pdf')),
           g, width = 14, height = 12)
    cat('  图已存到', output_path(paste0('fig1/', s, '_tsne_panels.pdf')), '\n')
}

# %%
main <- function() {
    # 允许命令行指定样本，便于单样本测试：
    #   Rscript projectmd/fig01_tsne_cnv_cytotrace_seurat.R GTE001
    args <- commandArgs(trailingOnly = TRUE)
    samples <- if (length(args) > 0) args else SAMPLES
    for (s in samples) {
        obj <- process_one(s)
        plots_for_sample(obj, s)
    }
    cat('\n✓ 全部完成\n')
}

# %%
main()
