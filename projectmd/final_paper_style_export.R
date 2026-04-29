#!/usr/bin/env Rscript
# final_paper_style_export.R — 把已有结果重新出图，逼近论文 Wu 2022 风格
# 输出：output/fig_paper_style/  （不动原 fig1/3/5）
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', script_args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
    script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
    setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')
options(Seurat.object.assay.version='v3')
suppressPackageStartupMessages({
    library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
    library(viridis); library(RColorBrewer); library(scales)
})
out <- output_path('fig_paper_style'); dir.create(out, showWarnings=FALSE, recursive=TRUE)

# ============================================================
# Fig 1A/B/D/E — GTE009 tSNE 4 panels
# A: CNV score (red gradient)
# B: Undifferentiated score (CytoTRACE, teal gradient)
# D: CNV subclone (3 colors)
# E: Cell type (10 colors)
# ============================================================
s <- readRDS(output_path('fig1/GTE009_seurat.rds'))
emb <- as.data.frame(Embeddings(s, 'tsne'))
colnames(emb) <- c('tsne_1','tsne_2')
emb <- cbind(emb, s@meta.data[, c('CNV_level','CytoTRACE','CNV_Cluster','Final_cluster')])
emb$CNV_Cluster <- as.character(emb$CNV_Cluster)
emb$CNV_Cluster[is.na(emb$CNV_Cluster)] <- 'NM'
emb$CNV_Cluster <- factor(emb$CNV_Cluster,
                          levels=c('Subclone_1','Subclone_2','NM'),
                          labels=c('1','2','NM'))

base_theme <- theme_void(base_size=11) +
    theme(plot.title = element_text(hjust=0.5, face='bold', size=12),
          legend.position = 'right',
          legend.key.size = unit(0.4,'cm'))

# Fig 1A: CNV score
p1A <- ggplot(emb, aes(tsne_1, tsne_2, color=CNV_level)) +
    geom_point(size=0.3, alpha=0.7) +
    scale_color_gradient(low='#FFE5E5', high='#B30000', na.value='grey85',
                         name='CNV\nscore', limits=c(0, 400), oob=squish) +
    ggtitle('CNV score') + base_theme

# Fig 1B: CytoTRACE
p1B <- ggplot(emb, aes(tsne_1, tsne_2, color=CytoTRACE)) +
    geom_point(size=0.3, alpha=0.7) +
    scale_color_gradient(low='grey90', high='#1B9E77', name='Undiff.\nscore') +
    ggtitle('Undifferentiated score') + base_theme

# Fig 1D: CNV subclone
p1D <- ggplot(emb, aes(tsne_1, tsne_2, color=CNV_Cluster)) +
    geom_point(size=0.3, alpha=0.7) +
    scale_color_manual(values=c('1'='#E64B35','2'='#4DBBD5','NM'='grey80'),
                       name='CNV\nsubclone') +
    ggtitle('CNV subclone') + base_theme +
    guides(color = guide_legend(override.aes = list(size=3)))

# Fig 1E: Cell type
ct_pal <- c('NSC-like'='#E64B35','RGC-like'='#9C27B0','Epe-like'='#2E7D32',
            'Neu-like'='#1976D2','Oli-like'='#FF9800','OPC-like'='#FBC02D',
            'Ast-like'='#00897B','OPC'='#FF6F61','Mic'='#7E57C2')
p1E <- ggplot(emb, aes(tsne_1, tsne_2, color=Final_cluster)) +
    geom_point(size=0.3, alpha=0.7) +
    scale_color_manual(values=ct_pal, name='Cell type') +
    ggtitle('GTE009 cell types') + base_theme +
    guides(color = guide_legend(override.aes = list(size=3), ncol=1))

ggsave(file.path(out,'Fig1A_CNV_score.pdf'), p1A, width=4.5, height=4)
ggsave(file.path(out,'Fig1B_Undiff_score.pdf'), p1B, width=4.5, height=4)
ggsave(file.path(out,'Fig1D_CNV_subclone.pdf'), p1D, width=4.5, height=4)
ggsave(file.path(out,'Fig1E_celltype.pdf'), p1E, width=5, height=4)
ggsave(file.path(out,'Fig1_panels_ABDE.pdf'),
       (p1A | p1B) / (p1D | p1E), width=10, height=8)
cat('✓ Fig 1 A/B/D/E\n')

# ============================================================
# Fig 1G — stacked bar of cell type per subclone (论文配色)
# ============================================================
df <- emb %>%
    filter(CNV_Cluster %in% c('1','2')) %>%
    count(CNV_Cluster, Final_cluster) %>%
    group_by(CNV_Cluster) %>% mutate(prop = n / sum(n))

p1G <- ggplot(df, aes(x=CNV_Cluster, y=prop, fill=Final_cluster)) +
    geom_bar(stat='identity', width=0.6, color='white', size=0.2) +
    scale_fill_manual(values=ct_pal, name='Cell type') +
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
    labs(x='Subclone', y='', title='Cell type composition by subclone') +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(hjust=0.5, face='bold'),
          legend.key.size=unit(0.4,'cm'))
ggsave(file.path(out,'Fig1G_stacked_bar.pdf'), p1G, width=5, height=5)
cat('✓ Fig 1G stacked bar\n')

# ============================================================
# Fig 3C — tSNE colored by trajectory score (red→green gradient)
# ============================================================
sig <- readRDS(output_path('fig3/trajectory_signatures.rds'))
expr <- GetAssayData(s, layer='data')
common <- intersect(rownames(expr), c(sig$undiff, sig$diff))
undiff_g <- intersect(sig$undiff, rownames(expr))
diff_g   <- intersect(sig$diff,   rownames(expr))
emb$traj_score <- as.numeric(Matrix::colMeans(expr[undiff_g, , drop=FALSE]) -
                             Matrix::colMeans(expr[diff_g, , drop=FALSE]))
# 标准化到 -1..1
emb$traj_score <- 2 * (emb$traj_score - min(emb$traj_score)) /
                  (max(emb$traj_score) - min(emb$traj_score)) - 1

p3C <- ggplot(emb, aes(tsne_1, tsne_2, color=traj_score)) +
    geom_point(size=0.3, alpha=0.8) +
    scale_color_gradient2(low='#1B9E77', mid='grey90', high='#E64B35',
                          midpoint=0, name='Trajectory\nscore',
                          limits=c(-1,1)) +
    ggtitle('Trajectory score (GTE009)') + base_theme
ggsave(file.path(out,'Fig3C_trajectory_score_tsne.pdf'), p3C, width=4.5, height=4)
cat('✓ Fig 3C trajectory tSNE\n')

# Violin per cell type
df_v <- emb %>% filter(!is.na(Final_cluster))
p3C_v <- ggplot(df_v, aes(x=Final_cluster, y=traj_score, fill=Final_cluster)) +
    geom_violin(scale='width', trim=TRUE, color=NA, alpha=0.8) +
    geom_boxplot(width=0.1, fill='white', outlier.shape=NA) +
    scale_fill_manual(values=ct_pal, guide='none') +
    labs(x='', y='Trajectory score', title='Trajectory score per cell type') +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(hjust=0.5, face='bold'),
          axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(out,'Fig3C_trajectory_violin.pdf'), p3C_v, width=6, height=4.5)
cat('✓ Fig 3C trajectory violin\n')

# ============================================================
# Fig 5A — donut chart of cell type counts (4 EPN samples 合计)
# ============================================================
# 重新从 4 个 seurat 拼总细胞型数
all_ct <- list()
for (s_id in c('GTE001','GTE002','GTE009','GTE012')) {
    f <- output_path(paste0('fig1/', s_id, '_seurat.rds'))
    if (!file.exists(f)) next
    so <- readRDS(f)
    # 优先 Final_cluster；没有就 Brief_cluster；再不然 Step1_cell_type
    col <- intersect(c('Final_cluster','Brief_cluster','Step1_cell_type'),
                     colnames(so@meta.data))[1]
    if (is.na(col)) next
    all_ct[[s_id]] <- table(so@meta.data[[col]])
}
ct_df <- bind_rows(lapply(all_ct, function(x) as.data.frame(x)), .id='sample')
colnames(ct_df) <- c('sample','celltype','n')
ct_total <- ct_df %>% group_by(celltype) %>% summarise(n=sum(n)) %>%
    arrange(desc(n)) %>% mutate(prop = n / sum(n))

# 论文 5A 是 donut；用 coord_polar
ct_total$ymax <- cumsum(ct_total$prop)
ct_total$ymin <- c(0, head(ct_total$ymax, -1))
ct_total$label <- paste0(ct_total$celltype, '\n', ct_total$n)

# 颜色：尽量沿用 ct_pal，未定义类别补色
pal_extra <- setNames(brewer.pal(8,'Set3'),
                      setdiff(as.character(ct_total$celltype), names(ct_pal)))
all_pal <- c(ct_pal, pal_extra)

p5A <- ggplot(ct_total, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.5, fill=celltype)) +
    geom_rect(color='white') +
    coord_polar(theta='y') + xlim(c(0,4)) +
    scale_fill_manual(values=all_pal, name='Cell type') +
    annotate('text', x=0, y=0, label=paste0('Total\n',sum(ct_total$n)),
             fontface='bold', size=4) +
    theme_void(base_size=11) +
    ggtitle('Cell numbers — 4 EPN samples') +
    theme(plot.title=element_text(hjust=0.5, face='bold'))
ggsave(file.path(out,'Fig5A_donut.pdf'), p5A, width=6, height=5)
cat('✓ Fig 5A donut\n')

# ============================================================
# Fig 1F — bulk DEG heatmap (top up + top dn 各 30)
# ============================================================
deg <- read.csv(output_path('fig1_bulk/bulk_DEG_relapse_vs_primary.csv'))
deg <- deg[!is.na(deg$adj.P.Val), ]
top_up <- head(deg[order(-deg$logFC), 'gene'], 30)
top_dn <- head(deg[order( deg$logFC), 'gene'], 30)
top_g  <- c(top_up, top_dn)

# 重新加载表达 + group 信息
suppressPackageStartupMessages({ library(GEOquery); library(hgu133plus2.db) })
gse <- getGEO(filename=raw_data_path('bulk/GSE64415_series_matrix.txt.gz'),
              getGPL=FALSE)
expr <- exprs(gse); pd <- pData(gse)
chars <- pd[, grep('characteristics', colnames(pd))]
extract <- function(key) apply(chars, 1, function(r) {
    x <- grep(key, r, value=TRUE, ignore.case=TRUE)
    if (length(x)) sub('.*: *','',x[1]) else NA
})
label   <- extract('primary/relapse'); subgrp  <- extract('molecular subgroup')
keep <- !is.na(label) & label %in% c('primary','relapse')
expr <- expr[, keep]; group <- factor(label[keep], levels=c('primary','relapse'))
subgrp <- subgrp[keep]
if (max(expr,na.rm=TRUE)>50) expr <- log2(expr+1)
expr <- expr[rowSums(is.na(expr))==0,]
map <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(expr),
                             columns='SYMBOL', keytype='PROBEID')
map <- map[!is.na(map$SYMBOL) & !duplicated(map$PROBEID),]
expr <- expr[map$PROBEID,]; expr <- limma::avereps(expr, ID=map$SYMBOL)

mat <- expr[intersect(top_g, rownames(expr)), ]
mat_z <- t(scale(t(mat)))
mat_z[mat_z >  2] <-  2; mat_z[mat_z < -2] <- -2
ord <- order(group, subgrp)
suppressPackageStartupMessages(library(pheatmap))
ann_col <- data.frame(Group=group[ord], Subgroup=subgrp[ord], row.names=colnames(mat_z)[ord])
ann_colors <- list(Group=c(primary='#4DBBD5', relapse='#E64B35'))
pdf(file.path(out,'Fig1F_bulk_heatmap.pdf'), width=10, height=8)
pheatmap(mat_z[, ord], cluster_cols=FALSE, cluster_rows=TRUE,
         annotation_col=ann_col, annotation_colors=ann_colors,
         show_colnames=FALSE, fontsize_row=6,
         color=colorRampPalette(c('#4575B4','white','#D73027'))(100),
         main='Bulk DEG: relapse vs primary (top 30↑ + 30↓)')
dev.off()
cat('✓ Fig 1F bulk heatmap\n')

cat('\n✅ 所有 paper-style 图完成 →', out, '\n')
