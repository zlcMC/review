#!/usr/bin/env Rscript
# Best-effort Figure 2C reproduction using Slingshot.
# Rationale: Monocle2/DDRTree is not stable in the current R+dplyr+igraph stack.
# This script keeps the biological target of Fig 2C: GTE009 malignant-cell trajectory.

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', script_args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
    script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
    setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')
options(Seurat.object.assay.version = 'v3')

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(slingshot)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    library(scales)
})

fig2_dir <- output_path('fig2')
paper_dir <- output_path('fig_paper_style')
dir.create(fig2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(output_path('fig1/GTE009_seurat.rds'))
mal <- subset(obj, subset = CNV_Cluster %in% c('Subclone_1', 'Subclone_2'))
cat('GTE009 malignant cells:', ncol(mal), '\n')

sce <- as.SingleCellExperiment(mal, assay = 'RNA')
cluster_labels <- factor(mal$Final_cluster)

# PCA is used for the analytical pseudotime because it is less distorted than tSNE.
sds_pca <- slingshot(
    sce,
    reducedDim = 'PCA',
    clusterLabels = cluster_labels,
    start.clus = 'NSC-like'
)
saveRDS(sds_pca, output_path('fig2/GTE009_fig2c_best_slingshot_pca.rds'))

lineages <- slingLineages(sds_pca)
lineage_end <- vapply(lineages, function(x) tail(x, 1), character(1))
main_lineage <- which(lineage_end == 'Epe-like')[1]
if (is.na(main_lineage)) main_lineage <- 1
cat('Main lineage:', main_lineage, paste(lineages[[main_lineage]], collapse = ' -> '), '\n')

pseudotime_matrix <- slingPseudotime(sds_pca)
main_pseudotime <- pseudotime_matrix[, main_lineage]
main_pseudotime_scaled <- rescale(main_pseudotime, to = c(0, 1),
                                  from = range(main_pseudotime, na.rm = TRUE))

# Fit a separate tSNE-space Slingshot object only for drawing smooth curves on the
# same coordinates used in Fig 1. The PCA pseudotime above remains the reported score.
sds_tsne <- slingshot(
    sce,
    reducedDim = 'TSNE',
    clusterLabels = cluster_labels,
    start.clus = 'NSC-like'
)
saveRDS(sds_tsne, output_path('fig2/GTE009_fig2c_best_slingshot_tsne.rds'))

embedding <- as.data.frame(Embeddings(mal, 'tsne'))
colnames(embedding) <- c('tsne_1', 'tsne_2')
embedding$cell <- rownames(embedding)

plot_data <- cbind(
    embedding,
    mal@meta.data[embedding$cell, c('Final_cluster', 'CNV_Cluster', 'CytoTRACE', 'S.Score', 'G2M.Score')]
)
plot_data$pseudotime <- main_pseudotime_scaled[embedding$cell]

scores_path <- output_path('fig3/GTE009_scores.rds')
if (file.exists(scores_path)) {
    scores <- readRDS(scores_path)
    common_cells <- intersect(rownames(scores), plot_data$cell)
    plot_data$trajectory_score <- NA_real_
    plot_data$trajectory_score[match(common_cells, plot_data$cell)] <- scores[common_cells, 'trajectory_score']
}

make_curve_df <- function(sds, main_lineage) {
    bind_rows(lapply(seq_along(slingCurves(sds)), function(lineage_id) {
        curve <- slingCurves(sds)[[lineage_id]]
        as.data.frame(curve$s[, 1:2, drop = FALSE]) %>%
            setNames(c('tsne_1', 'tsne_2')) %>%
            mutate(
                lineage = paste0('Lineage ', lineage_id),
                is_main = lineage_id == main_lineage,
                order = row_number()
            )
    }))
}
curve_data <- make_curve_df(sds_tsne, main_lineage)
all_curves <- curve_data %>% filter(!is_main)
main_curve <- curve_data %>% filter(is_main)

celltype_palette <- c(
    'NSC-like' = '#D73027',
    'RGC-like' = '#984EA3',
    'Ast-like' = '#00897B',
    'Epe-like' = '#2E7D32',
    'OPC-like' = '#FBC02D',
    'Oli-like' = '#FF9800',
    'Neu-like' = '#1976D2'
)
subclone_palette <- c('Subclone_1' = '#E64B35', 'Subclone_2' = '#4DBBD5')

base_theme <- theme_void(base_size = 11) +
    theme(
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 12),
        legend.position = 'right',
        legend.key.size = unit(0.35, 'cm')
    )

lineage_layers <- list(
    geom_path(
        data = all_curves,
        aes(tsne_1, tsne_2, group = lineage),
        inherit.aes = FALSE,
        color = 'grey55', linewidth = 0.35, alpha = 0.45
    ),
    geom_path(
        data = main_curve,
        aes(tsne_1, tsne_2),
        inherit.aes = FALSE,
        color = 'black', linewidth = 0.8,
        arrow = arrow(length = unit(0.08, 'inches'), type = 'closed')
    )
)

label_data <- plot_data %>%
    group_by(Final_cluster) %>%
    summarise(tsne_1 = median(tsne_1), tsne_2 = median(tsne_2), .groups = 'drop')

p_pseudotime <- ggplot(plot_data, aes(tsne_1, tsne_2, color = pseudotime)) +
    geom_point(size = 0.28, alpha = 0.85, na.rm = FALSE) +
    lineage_layers +
    scale_color_gradient(low = '#2C7BB6', high = '#D7191C', na.value = 'grey88',
                         limits = c(0, 1), name = 'Pseudotime') +
    ggtitle('Pseudotime') + base_theme

p_state <- ggplot(plot_data, aes(tsne_1, tsne_2, color = Final_cluster)) +
    geom_point(size = 0.28, alpha = 0.85) +
    lineage_layers +
    geom_text(data = label_data, aes(x = tsne_1, y = tsne_2, label = Final_cluster), inherit.aes = FALSE,
              size = 2.7, fontface = 'bold', color = 'black') +
    scale_color_manual(values = celltype_palette, name = 'Cell state', na.value = 'grey85') +
    ggtitle('Cell state') + base_theme +
    guides(color = guide_legend(override.aes = list(size = 2.5)))

p_subclone <- ggplot(plot_data, aes(tsne_1, tsne_2, color = CNV_Cluster)) +
    geom_point(size = 0.28, alpha = 0.85) +
    lineage_layers +
    scale_color_manual(values = subclone_palette, name = 'Subclone', na.value = 'grey85') +
    ggtitle('CNV subclone') + base_theme +
    guides(color = guide_legend(override.aes = list(size = 2.5)))

p_cytotrace <- ggplot(plot_data, aes(tsne_1, tsne_2, color = CytoTRACE)) +
    geom_point(size = 0.28, alpha = 0.85) +
    lineage_layers +
    scale_color_gradient(low = 'grey90', high = '#1B9E77', name = 'CytoTRACE') +
    ggtitle('CytoTRACE') + base_theme

panel <- (p_pseudotime | p_state) / (p_subclone | p_cytotrace)

ggsave(output_path('fig2/Fig2C_GTE009_best_slingshot.pdf'), panel, width = 10, height = 8)
ggsave(file.path(paper_dir, 'Fig2C_GTE009_best_slingshot.pdf'), panel, width = 10, height = 8)
ggsave(output_path('fig2/Fig2C_GTE009_best_slingshot.png'), panel, width = 10, height = 8, dpi = 300)

single_panel <- p_pseudotime + ggtitle('GTE009 malignant trajectory')
ggsave(file.path(paper_dir, 'Fig2C_GTE009_best_slingshot_pseudotime_only.pdf'),
       single_panel, width = 5.2, height = 4.6)

summary_by_state <- plot_data %>%
    group_by(Final_cluster) %>%
    summarise(
        n = n(),
        pseudotime_mean = mean(pseudotime, na.rm = TRUE),
        pseudotime_median = median(pseudotime, na.rm = TRUE),
        CytoTRACE_mean = mean(CytoTRACE, na.rm = TRUE),
        Subclone_1 = sum(CNV_Cluster == 'Subclone_1', na.rm = TRUE),
        Subclone_2 = sum(CNV_Cluster == 'Subclone_2', na.rm = TRUE),
        .groups = 'drop'
    )
write.csv(summary_by_state, output_path('fig2/Fig2C_best_summary_by_state.csv'), row.names = FALSE)

summary_by_subclone <- plot_data %>%
    group_by(CNV_Cluster) %>%
    summarise(
        n = n(),
        pseudotime_mean = mean(pseudotime, na.rm = TRUE),
        pseudotime_median = median(pseudotime, na.rm = TRUE),
        CytoTRACE_mean = mean(CytoTRACE, na.rm = TRUE),
        S_mean = mean(S.Score, na.rm = TRUE),
        G2M_mean = mean(G2M.Score, na.rm = TRUE),
        .groups = 'drop'
    )
write.csv(summary_by_subclone, output_path('fig2/Fig2C_best_summary_by_subclone.csv'), row.names = FALSE)

write.csv(
    plot_data[, intersect(c('cell', 'Final_cluster', 'CNV_Cluster', 'pseudotime',
                            'CytoTRACE', 'S.Score', 'G2M.Score', 'trajectory_score'),
                          colnames(plot_data))],
    output_path('fig2/Fig2C_best_cell_scores.csv'),
    row.names = FALSE
)

cat('OK: Fig2C best Slingshot outputs written to output/fig2 and output/fig_paper_style\n')