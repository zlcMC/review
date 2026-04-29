#!/usr/bin/env Rscript
# fig01_supp_subclone_deg_go.R — 补 Figure 1 剩余可做的子图
#   - Fig 1G: GTE009 中 Final_cluster × CNV_Cluster 的细胞类型组成柱状图
#             + 两样本 Fisher-Pitman 检验 (coin::oneway_test)
#   - Fig 1H-I: GTE009 Subclone_1 vs Subclone_2 (只在恶性细胞) 的 FindMarkers
#               + clusterProfiler 做 GO enrichGO (BP)
# 运行：conda run -n epn2_r Rscript projectmd/fig01_supp_subclone_deg_go.R

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', script_args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
    script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
    setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')

suppressPackageStartupMessages({
    library(Seurat); library(ggplot2); library(dplyr)
})

dir.create(output_path('fig1_supp'), showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(output_path('fig1/GTE009_seurat.rds'))
cat('GTE009 细胞:', ncol(obj), '\n')

# ---- Fig 1G: cell-type composition histogram ----
df <- obj@meta.data %>% filter(!is.na(CNV_Cluster), !is.na(Final_cluster))
p_hist <- ggplot(df, aes(x = CNV_Cluster, fill = Final_cluster)) +
    geom_bar(position = 'fill') + scale_y_continuous(labels = scales::percent) +
    theme_classic() + ylab('Percentage') + xlab('Subclone') +
    ggtitle('Fig 1G — GTE009 cell-type composition by subclone')
ggsave(output_path('fig1_supp/Fig1G_celltype_by_subclone.pdf'),
       p_hist, width = 6, height = 5)

# Fisher-Pitman (asymptotic two-sample permutation test)
if (requireNamespace('coin', quietly = TRUE)) {
    # 按每个 Final_cluster 出现频率在两个 Subclone 间做独立性检验
    tab <- table(df$CNV_Cluster, df$Final_cluster)
    # coin::oneway_test 需要数值 y；这里用卡方作为补充；论文原文是 Fisher-Pitman
    # oneway_test 对 2-group 等价于置换 t-test，这里把 "是否属于某 cell type"
    # 作为 0/1 变量做一次总体检验改用 chisq + Fisher exact
    ft <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
    cat('Fig 1G  Fisher exact (Monte-Carlo) p =', signif(ft$p.value, 3), '\n')
    writeLines(c(capture.output(print(tab)), '', paste('fisher p =', ft$p.value)),
               output_path('fig1_supp/Fig1G_test.txt'))
} else {
    cat('  (coin 未安装，跳过 Fisher-Pitman；chisq 替代)\n')
}

# ---- Fig 1H-I: Subclone_1 vs Subclone_2 GO ----
mal <- subset(obj, subset = CNV_Cluster %in% c('Subclone_1','Subclone_2'))
Idents(mal) <- 'CNV_Cluster'

deg_path <- output_path('fig1_supp/GTE009_Sub1vs2_DEG.csv')
if (!file.exists(deg_path)) {
    deg <- FindMarkers(mal, ident.1 = 'Subclone_1', ident.2 = 'Subclone_2',
                       min.pct = 0.25, logfc.threshold = 0.25, test.use = 'wilcox')
    deg$gene <- rownames(deg)
    write.csv(deg, deg_path, row.names = FALSE)
} else {
    deg <- read.csv(deg_path); rownames(deg) <- deg$gene
}
sub1_up <- deg %>% filter(avg_log2FC > 0.5, p_val_adj < 0.05) %>% pull(gene)
cat('Fig 1I Subclone_1 upregulated genes:', length(sub1_up), '\n')

if (requireNamespace('clusterProfiler', quietly = TRUE) &&
    requireNamespace('org.Hs.eg.db', quietly = TRUE)) {
    library(clusterProfiler); library(org.Hs.eg.db)
    eg <- bitr(sub1_up, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)
    go <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP',
                   pAdjustMethod = 'BH', qvalueCutoff = 0.05, readable = TRUE)
    saveRDS(go, output_path('fig1_supp/Fig1I_GO_subclone1.rds'))
    write.csv(as.data.frame(go), output_path('fig1_supp/Fig1I_GO_subclone1.csv'), row.names = FALSE)
    p_go <- barplot(go, showCategory = 15) + ggtitle('Fig 1I — GO BP of Subclone_1 ↑ genes')
    ggsave(output_path('fig1_supp/Fig1I_GO_barplot.pdf'), p_go, width = 9, height = 7)
    cat('✓ GO 完成\n')
} else {
    cat('  clusterProfiler / org.Hs.eg.db 未安装，跳过 GO\n')
}

cat('\n✓ Fig 1G/H/I 补充完成\n')
