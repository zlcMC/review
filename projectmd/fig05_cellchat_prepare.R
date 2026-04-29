#!/usr/bin/env Rscript
# fig05_cellchat_prepare.R (v2) — 内存友好版本
# 论文 Fig 5 + Supp 8: CellChat crosstalk on integrated 4 EPN samples
# 改动：放弃 Seurat 5 merge (会 dense 化导致 OOM)
#       改用手动 sparse cbind 拼接 4 样本的 RNA counts，
#       直接在 sparse 矩阵上做 LogNormalize，再喂给 CellChat。
#       全程保持 sparse，对 30K 细胞 × 20K 基因，内存 < 2GB。
# 运行：conda run -n epn2_r Rscript projectmd/fig05_cellchat_prepare.R

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', script_args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
    script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
    setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')

suppressPackageStartupMessages({
    library(Seurat); library(Matrix); library(CellChat)
    library(ggplot2); library(dplyr); library(patchwork)
})

SAMPLES <- c('GTE001','GTE002','GTE009','GTE012')
dir.create(output_path('fig5'), showWarnings = FALSE, recursive = TRUE)

# ---- 1) sparse cbind 拼接 ----
cache <- output_path('fig5/merged_sparse.rds')
if (file.exists(cache)) {
    cat('载入缓存\n'); merged <- readRDS(cache)
} else {
    counts_list <- list(); meta_list <- list()
    for (s in SAMPLES) {
        o <- readRDS(output_path(paste0('fig1/', s, '_seurat.rds')))
        cm <- GetAssayData(o, assay='RNA', layer='counts')
        colnames(cm) <- paste0(s, '_', colnames(cm))
        counts_list[[s]] <- cm
        m <- o@meta.data; rownames(m) <- paste0(s, '_', rownames(m))
        m$sample <- s
        meta_list[[s]] <- m
        rm(o); invisible(gc())
        cat(sprintf('  %s: %d 基因 × %d 细胞\n', s, nrow(cm), ncol(cm)))
    }
    common <- Reduce(intersect, lapply(counts_list, rownames))
    cat(sprintf('  公共基因: %d\n', length(common)))
    counts_list <- lapply(counts_list, function(x) x[common, ])
    all_counts <- do.call(cbind, counts_list)
    rm(counts_list); invisible(gc())
    common_meta_cols <- Reduce(intersect, lapply(meta_list, colnames))
    all_meta <- do.call(rbind, unname(lapply(meta_list, function(m) m[, common_meta_cols, drop=FALSE])))
    stopifnot(identical(colnames(all_counts), rownames(all_meta)))
    cat(sprintf('  合并后: %d 细胞 × %d 基因\n', ncol(all_counts), nrow(all_counts)))
    sf <- Matrix::colSums(all_counts)
    all_data <- all_counts %*% Matrix::Diagonal(x = 1e4 / sf)
    all_data@x <- log1p(all_data@x)
    dimnames(all_data) <- dimnames(all_counts)
    merged <- list(counts = all_counts, data = all_data, meta = all_meta)
    saveRDS(merged, cache)
}
cat(sprintf('合并: %d 细胞\n', ncol(merged$data)))

# ---- 2) 导出给 pySCENIC ----
scenic_dir <- output_path('fig6/scenic_input')
dir.create(scenic_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(expr = merged$data, meta = merged$meta),
        file.path(scenic_dir, 'integrated_logmat.rds'))
cat('已存 pySCENIC 输入\n')

# ---- 3) Fig 5A: cell-type 数 per sample ----
group_col <- if ('Final_cluster' %in% colnames(merged$meta)) 'Final_cluster' else 'Brief_cluster'
cnt_df <- as.data.frame(merged$meta)
cnt_df$.grp <- cnt_df[[group_col]]
cnt <- cnt_df %>% count(sample, .grp)
p_cnt <- ggplot(cnt, aes(x = sample, y = n, fill = .grp)) +
         geom_bar(stat = 'identity') + theme_classic() + ylab('# cells') +
         labs(fill = group_col,
              title = paste('Fig 5A - cell numbers per sample (group =', group_col, ')'))
ggsave(output_path('fig5/Fig5A_cell_numbers.pdf'), p_cnt, width = 8, height = 5)

# ---- 4) CellChat ----
cc_cache <- output_path('fig5/cellchat_merged.rds')
if (file.exists(cc_cache)) {
    cellchat <- readRDS(cc_cache)
} else {
    labels <- as.character(merged$meta[[group_col]])
    keep   <- !is.na(labels) & labels != ''
    data_use <- merged$data[, keep]
    meta_cc  <- data.frame(group = labels[keep],
                           sample = merged$meta$sample[keep],
                           samples = merged$meta$sample[keep],
                           row.names = colnames(data_use))
    cellchat <- createCellChat(object = data_use, meta = meta_cc, group.by = 'group')
    cellchat@DB <- CellChatDB.human
    cellchat <- subsetData(cellchat)
    future::plan('sequential')
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    saveRDS(cellchat, cc_cache)
}

cat('Fig 5A and CellChat object done. Run fig05_cellchat_downstream.R for final Fig 5/Supp Fig 8 figures.\n')
