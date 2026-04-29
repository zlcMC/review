#!/usr/bin/env Rscript
# v2 of final_paper_style_export.R — paper-style summary figures.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

options(Seurat.object.assay.version = 'v3')
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
  library(viridis); library(RColorBrewer); library(scales)
})
out <- v2_dir('fig_paper_style')

# ---- Fig 1 A/B/D/E: GTE009 tSNE ----
s <- readRDS(dep_file('fig1', 'GTE009_seurat.rds'))
emb <- as.data.frame(Embeddings(s, 'tsne'))
colnames(emb) <- c('tsne_1','tsne_2')
emb <- cbind(emb, s@meta.data[, c('CNV_level','CytoTRACE','CNV_Cluster','Final_cluster')])
emb$CNV_Cluster <- as.character(emb$CNV_Cluster)
emb$CNV_Cluster[is.na(emb$CNV_Cluster)] <- 'NM'
emb$CNV_Cluster <- factor(emb$CNV_Cluster,
                          levels = c('Subclone_1','Subclone_2','NM'),
                          labels = c('1','2','NM'))

base_theme <- theme_void(base_size = 11) +
  theme(plot.title      = element_text(hjust = 0.5, face = 'bold', size = 12),
        legend.position = 'right',
        legend.key.size = unit(0.4, 'cm'))

ct_pal <- c('NSC-like'='#E64B35','RGC-like'='#9C27B0','Epe-like'='#2E7D32',
            'Neu-like'='#1976D2','Oli-like'='#FF9800','OPC-like'='#FBC02D',
            'Ast-like'='#00897B','OPC'='#FF6F61','Mic'='#7E57C2')

tsne_panel <- function(color_aes, scale_layer, title, file, w = 4.5, h = 4) {
  p <- ggplot(emb, aes(tsne_1, tsne_2)) +
    geom_point(aes(color = !!color_aes), size = 0.3, alpha = 0.7) +
    scale_layer + ggtitle(title) + base_theme
  ggsave(file.path(out, file), p, width = w, height = h)
  p
}

p1A <- tsne_panel(quote(CNV_level),
                  scale_color_gradient(low = '#FFE5E5', high = '#B30000',
                                       na.value = 'grey85', name = 'CNV\nscore',
                                       limits = c(0, 400), oob = squish),
                  'CNV score', 'Fig1A_CNV_score.pdf')
p1B <- tsne_panel(quote(CytoTRACE),
                  scale_color_gradient(low = 'grey90', high = '#1B9E77',
                                       name = 'Undiff.\nscore'),
                  'Undifferentiated score', 'Fig1B_Undiff_score.pdf')
p1D <- tsne_panel(quote(CNV_Cluster),
                  list(scale_color_manual(
                         values = c('1' = '#E64B35', '2' = '#4DBBD5', 'NM' = 'grey80'),
                         name = 'CNV\nsubclone'),
                       guides(color = guide_legend(override.aes = list(size = 3)))),
                  'CNV subclone', 'Fig1D_CNV_subclone.pdf')
p1E <- tsne_panel(quote(Final_cluster),
                  list(scale_color_manual(values = ct_pal, name = 'Cell type'),
                       guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))),
                  'GTE009 cell types', 'Fig1E_celltype.pdf', w = 5)
ggsave(file.path(out, 'Fig1_panels_ABDE.pdf'),
       (p1A | p1B) / (p1D | p1E), width = 10, height = 8)
cat('✓ Fig 1 A/B/D/E\n')

# ---- Fig 1G: stacked bar ----
df <- emb %>% filter(CNV_Cluster %in% c('1','2')) %>%
  count(CNV_Cluster, Final_cluster) %>%
  group_by(CNV_Cluster) %>% mutate(prop = n / sum(n))
ggsave(file.path(out, 'Fig1G_stacked_bar.pdf'),
       ggplot(df, aes(x = CNV_Cluster, y = prop, fill = Final_cluster)) +
         geom_bar(stat = 'identity', width = 0.6, color = 'white', size = 0.2) +
         scale_fill_manual(values = ct_pal, name = 'Cell type') +
         scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
         labs(x = 'Subclone', y = '', title = 'Cell type composition by subclone') +
         theme_classic(12) +
         theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
               legend.key.size = unit(0.4, 'cm')),
       width = 5, height = 5)
cat('✓ Fig 1G stacked bar\n')

# ---- Fig 3C: trajectory score tSNE + violin ----
sig <- readRDS(dep_file('fig3', 'trajectory_signatures.rds'))
expr <- GetAssayData(s, layer = 'data')
undiff_g <- intersect(sig$undiff, rownames(expr))
diff_g   <- intersect(sig$diff,   rownames(expr))
emb$traj_score <- as.numeric(Matrix::colMeans(expr[undiff_g, , drop = FALSE]) -
                             Matrix::colMeans(expr[diff_g,   , drop = FALSE]))
emb$traj_score <- 2 * (emb$traj_score - min(emb$traj_score)) /
                  (max(emb$traj_score) - min(emb$traj_score)) - 1

ggsave(file.path(out, 'Fig3C_trajectory_score_tsne.pdf'),
       ggplot(emb, aes(tsne_1, tsne_2, color = traj_score)) +
         geom_point(size = 0.3, alpha = 0.8) +
         scale_color_gradient2(low = '#1B9E77', mid = 'grey90', high = '#E64B35',
                               midpoint = 0, name = 'Trajectory\nscore',
                               limits = c(-1, 1)) +
         ggtitle('Trajectory score (GTE009)') + base_theme,
       width = 4.5, height = 4)

ggsave(file.path(out, 'Fig3C_trajectory_violin.pdf'),
       ggplot(emb %>% filter(!is.na(Final_cluster)),
              aes(x = Final_cluster, y = traj_score, fill = Final_cluster)) +
         geom_violin(scale = 'width', trim = TRUE, color = NA, alpha = 0.8) +
         geom_boxplot(width = 0.1, fill = 'white', outlier.shape = NA) +
         scale_fill_manual(values = ct_pal, guide = 'none') +
         labs(x = '', y = 'Trajectory score',
              title = 'Trajectory score per cell type') +
         theme_classic(11) +
         theme(plot.title  = element_text(hjust = 0.5, face = 'bold'),
               axis.text.x = element_text(angle = 45, hjust = 1)),
       width = 6, height = 4.5)
cat('✓ Fig 3C trajectory tSNE/violin\n')

# ---- Fig 5A donut: combine 4 EPN samples ----
all_ct <- list()
for (s_id in c('GTE001','GTE002','GTE009','GTE012')) {
  f <- dep_file('fig1', paste0(s_id, '_seurat.rds'))
  if (!file.exists(f)) next
  so <- readRDS(f)
  col <- intersect(c('Final_cluster','Brief_cluster','Step1_cell_type'),
                   colnames(so@meta.data))[1]
  if (is.na(col)) next
  all_ct[[s_id]] <- table(so@meta.data[[col]])
}
ct_df <- bind_rows(lapply(all_ct, function(x) as.data.frame(x)), .id = 'sample')
colnames(ct_df) <- c('sample','celltype','n')
ct_total <- ct_df %>% group_by(celltype) %>% summarise(n = sum(n)) %>%
  arrange(desc(n)) %>% mutate(prop = n / sum(n))
ct_total$ymax <- cumsum(ct_total$prop)
ct_total$ymin <- c(0, head(ct_total$ymax, -1))
pal_extra <- setNames(brewer.pal(8, 'Set3'),
                      setdiff(as.character(ct_total$celltype), names(ct_pal)))
all_pal <- c(ct_pal, pal_extra)
ggsave(file.path(out, 'Fig5A_donut.pdf'),
       ggplot(ct_total, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5,
                            fill = celltype)) +
         geom_rect(color = 'white') +
         coord_polar(theta = 'y') + xlim(c(0, 4)) +
         scale_fill_manual(values = all_pal, name = 'Cell type') +
         annotate('text', x = 0, y = 0,
                  label = paste0('Total\n', sum(ct_total$n)),
                  fontface = 'bold', size = 4) +
         theme_void(11) +
         ggtitle('Cell numbers — 4 EPN samples') +
         theme(plot.title = element_text(hjust = 0.5, face = 'bold')),
       width = 6, height = 5)
cat('✓ Fig 5A donut\n')

# ---- Fig 1F bulk DEG heatmap (uses load_gse_expression) ----
deg <- read.csv(dep_file('fig1_bulk', 'bulk_DEG_relapse_vs_primary.csv'))
deg <- deg[!is.na(deg$adj.P.Val), ]
top_g <- c(head(deg[order(-deg$logFC), 'gene'], 30),
           head(deg[order( deg$logFC), 'gene'], 30))

bulk <- load_gse_expression(
  matrix_path          = raw_data_path('bulk/GSE64415_series_matrix.txt.gz'),
  extra_columns        = list(label = 'primary/relapse', subgroup = 'molecular subgroup'),
  sample_filter_key    = 'label',
  sample_filter_values = c('primary','relapse')
)
expr  <- bulk$expr
group <- factor(bulk$meta$label, levels = c('primary','relapse'))

mat   <- expr[intersect(top_g, rownames(expr)), ]
mat_z <- t(scale(t(mat)))
mat_z[mat_z >  2] <-  2; mat_z[mat_z < -2] <- -2
ord   <- order(group, bulk$meta$subgroup)
suppressPackageStartupMessages(library(pheatmap))
ann_col    <- data.frame(Group = group[ord], Subgroup = bulk$meta$subgroup[ord],
                         row.names = colnames(mat_z)[ord])
ann_colors <- list(Group = c(primary = '#4DBBD5', relapse = '#E64B35'))
pdf(file.path(out, 'Fig1F_bulk_heatmap.pdf'), width = 10, height = 8)
pheatmap(mat_z[, ord], cluster_cols = FALSE, cluster_rows = TRUE,
         annotation_col = ann_col, annotation_colors = ann_colors,
         show_colnames = FALSE, fontsize_row = 6,
         color = colorRampPalette(c('#4575B4','white','#D73027'))(100),
         main = 'Bulk DEG: relapse vs primary (top 30↑ + 30↓)')
dev.off()
cat('✓ Fig 1F bulk heatmap\n')

cat('\n✅ paper-style 图完成 →', out, '\n')
