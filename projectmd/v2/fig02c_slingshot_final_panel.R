#!/usr/bin/env Rscript
# fig02c_slingshot_final_panel.R (v2) — 论文用 Fig2C：GTE009 恶性细胞 Slingshot 轨迹
# PCA 上跑 pseudotime（更稳），tSNE 上画曲线（与 Fig 1 一致）。
# 运行：conda run -n epn2_r Rscript projectmd/v2/fig02c_slingshot_final_panel.R

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()
options(Seurat.object.assay.version = 'v3')

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(slingshot)
  library(ggplot2); library(dplyr); library(patchwork); library(scales)
})

paper_dir  <- v2_dir('fig_paper_style')
fig2_dir   <- v2_dir('fig2')

obj <- readRDS(dep_file('fig1', 'GTE009_seurat.rds'))
mal <- subset(obj, subset = CNV_Cluster %in% c('Subclone_1', 'Subclone_2'))
cat('GTE009 malignant cells:', ncol(mal), '\n')

sce <- as.SingleCellExperiment(mal, assay = 'RNA')
clu <- factor(mal$Final_cluster)

sds_pca <- slingshot(sce, reducedDim = 'PCA', clusterLabels = clu, start.clus = 'NSC-like')
saveRDS(sds_pca, v2_file('fig2', 'GTE009_fig2c_best_slingshot_pca.rds'))

lineages    <- slingLineages(sds_pca)
end_state   <- vapply(lineages, function(x) tail(x, 1), character(1))
main_idx    <- which(end_state == 'Epe-like')[1]
if (is.na(main_idx)) main_idx <- 1
cat('Main lineage:', main_idx, paste(lineages[[main_idx]], collapse = ' -> '), '\n')

pt <- slingPseudotime(sds_pca)[, main_idx]
pt_scaled <- rescale(pt, to = c(0, 1), from = range(pt, na.rm = TRUE))

sds_tsne <- slingshot(sce, reducedDim = 'TSNE', clusterLabels = clu, start.clus = 'NSC-like')
saveRDS(sds_tsne, v2_file('fig2', 'GTE009_fig2c_best_slingshot_tsne.rds'))

emb <- as.data.frame(Embeddings(mal, 'tsne'))
colnames(emb) <- c('tsne_1', 'tsne_2'); emb$cell <- rownames(emb)

plot_data <- cbind(emb, mal@meta.data[emb$cell,
  c('Final_cluster', 'CNV_Cluster', 'CytoTRACE', 'S.Score', 'G2M.Score')])
plot_data$pseudotime <- pt_scaled[emb$cell]

scores_path <- dep_file('fig3', 'GTE009_scores.rds')
if (file.exists(scores_path)) {
  sc <- readRDS(scores_path)
  common <- intersect(rownames(sc), plot_data$cell)
  plot_data$trajectory_score <- NA_real_
  plot_data$trajectory_score[match(common, plot_data$cell)] <- sc[common, 'trajectory_score']
}

curve_df <- function(sds, main_idx) {
  bind_rows(lapply(seq_along(slingCurves(sds)), function(i) {
    s <- slingCurves(sds)[[i]]$s[, 1:2, drop = FALSE]
    as.data.frame(s) %>% setNames(c('tsne_1', 'tsne_2')) %>%
      mutate(lineage = paste0('Lineage ', i), is_main = i == main_idx, order = row_number())
  }))
}
curves <- curve_df(sds_tsne, main_idx)
others <- curves %>% filter(!is_main)
main_c <- curves %>% filter(is_main)

celltype_palette <- c('NSC-like'='#D73027','RGC-like'='#984EA3','Ast-like'='#00897B',
  'Epe-like'='#2E7D32','OPC-like'='#FBC02D','Oli-like'='#FF9800','Neu-like'='#1976D2')
subclone_palette <- c('Subclone_1' = '#E64B35', 'Subclone_2' = '#4DBBD5')

base_theme <- theme_void(base_size = 11) + theme(
  plot.title = element_text(hjust = 0.5, face = 'bold', size = 12),
  legend.position = 'right', legend.key.size = unit(0.35, 'cm'))

lineage_layers <- list(
  geom_path(data = others, aes(tsne_1, tsne_2, group = lineage),
            inherit.aes = FALSE, color = 'grey55', linewidth = 0.35, alpha = 0.45),
  geom_path(data = main_c, aes(tsne_1, tsne_2), inherit.aes = FALSE,
            color = 'black', linewidth = 0.8,
            arrow = arrow(length = unit(0.08, 'inches'), type = 'closed')))

label_data <- plot_data %>% group_by(Final_cluster) %>%
  summarise(tsne_1 = median(tsne_1), tsne_2 = median(tsne_2), .groups = 'drop')

p_pt <- ggplot(plot_data, aes(tsne_1, tsne_2, color = pseudotime)) +
  geom_point(size = 0.28, alpha = 0.85) + lineage_layers +
  scale_color_gradient(low = '#2C7BB6', high = '#D7191C', na.value = 'grey88',
                       limits = c(0, 1), name = 'Pseudotime') +
  ggtitle('Pseudotime') + base_theme

p_state <- ggplot(plot_data, aes(tsne_1, tsne_2, color = Final_cluster)) +
  geom_point(size = 0.28, alpha = 0.85) + lineage_layers +
  geom_text(data = label_data, aes(tsne_1, tsne_2, label = Final_cluster),
            inherit.aes = FALSE, size = 2.7, fontface = 'bold', color = 'black') +
  scale_color_manual(values = celltype_palette, name = 'Cell state', na.value = 'grey85') +
  ggtitle('Cell state') + base_theme +
  guides(color = guide_legend(override.aes = list(size = 2.5)))

p_sub <- ggplot(plot_data, aes(tsne_1, tsne_2, color = CNV_Cluster)) +
  geom_point(size = 0.28, alpha = 0.85) + lineage_layers +
  scale_color_manual(values = subclone_palette, name = 'Subclone', na.value = 'grey85') +
  ggtitle('CNV subclone') + base_theme +
  guides(color = guide_legend(override.aes = list(size = 2.5)))

p_cyto <- ggplot(plot_data, aes(tsne_1, tsne_2, color = CytoTRACE)) +
  geom_point(size = 0.28, alpha = 0.85) + lineage_layers +
  scale_color_gradient(low = 'grey90', high = '#1B9E77', name = 'CytoTRACE') +
  ggtitle('CytoTRACE') + base_theme

panel <- (p_pt | p_state) / (p_sub | p_cyto)
ggsave(v2_file('fig2', 'Fig2C_GTE009_best_slingshot.pdf'), panel, width = 10, height = 8)
ggsave(file.path(paper_dir, 'Fig2C_GTE009_best_slingshot.pdf'), panel, width = 10, height = 8)
ggsave(v2_file('fig2', 'Fig2C_GTE009_best_slingshot.png'), panel, width = 10, height = 8, dpi = 300)
ggsave(file.path(paper_dir, 'Fig2C_GTE009_best_slingshot_pseudotime_only.pdf'),
       p_pt + ggtitle('GTE009 malignant trajectory'), width = 5.2, height = 4.6)

summary_state <- plot_data %>% group_by(Final_cluster) %>%
  summarise(n = n(),
            pseudotime_mean = mean(pseudotime, na.rm = TRUE),
            pseudotime_median = median(pseudotime, na.rm = TRUE),
            CytoTRACE_mean = mean(CytoTRACE, na.rm = TRUE),
            Subclone_1 = sum(CNV_Cluster == 'Subclone_1', na.rm = TRUE),
            Subclone_2 = sum(CNV_Cluster == 'Subclone_2', na.rm = TRUE), .groups = 'drop')
write.csv(summary_state, v2_file('fig2', 'Fig2C_best_summary_by_state.csv'), row.names = FALSE)

summary_sub <- plot_data %>% group_by(CNV_Cluster) %>%
  summarise(n = n(),
            pseudotime_mean = mean(pseudotime, na.rm = TRUE),
            pseudotime_median = median(pseudotime, na.rm = TRUE),
            CytoTRACE_mean = mean(CytoTRACE, na.rm = TRUE),
            S_mean = mean(S.Score, na.rm = TRUE),
            G2M_mean = mean(G2M.Score, na.rm = TRUE), .groups = 'drop')
write.csv(summary_sub, v2_file('fig2', 'Fig2C_best_summary_by_subclone.csv'), row.names = FALSE)

write.csv(
  plot_data[, intersect(c('cell','Final_cluster','CNV_Cluster','pseudotime',
                          'CytoTRACE','S.Score','G2M.Score','trajectory_score'),
                        colnames(plot_data))],
  v2_file('fig2', 'Fig2C_best_cell_scores.csv'), row.names = FALSE)

cat('OK: Fig2C best Slingshot outputs written under output/v2/\n')
