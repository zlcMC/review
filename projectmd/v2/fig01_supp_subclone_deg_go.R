#!/usr/bin/env Rscript
# v2 of fig01_supp_subclone_deg_go.R
# Fig 1G/H/I: subclone composition + Subclone_1 vs Subclone_2 DEG/GO.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr)
})

obj <- readRDS(dep_file('fig1', 'GTE009_seurat.rds'))
cat('GTE009 cells:', ncol(obj), '\n')

# ---- Fig 1G: cell-type composition by subclone ----
df <- obj@meta.data %>% filter(!is.na(CNV_Cluster), !is.na(Final_cluster))
ggsave(
  v2_file('fig1_supp', 'Fig1G_celltype_by_subclone.pdf'),
  ggplot(df, aes(x = CNV_Cluster, fill = Final_cluster)) +
    geom_bar(position = 'fill') +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() + ylab('Percentage') + xlab('Subclone') +
    ggtitle('Fig 1G — GTE009 cell-type composition by subclone'),
  width = 6, height = 5
)

tab <- table(df$CNV_Cluster, df$Final_cluster)
ft  <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
cat('Fig 1G  Fisher exact (Monte-Carlo) p =', signif(ft$p.value, 3), '\n')
writeLines(c(capture.output(print(tab)), '', paste('fisher p =', ft$p.value)),
           v2_file('fig1_supp', 'Fig1G_test.txt'))

# ---- Fig 1H-I: Subclone_1 vs Subclone_2 DEG + GO ----
mal <- subset(obj, subset = CNV_Cluster %in% c('Subclone_1', 'Subclone_2'))
Idents(mal) <- 'CNV_Cluster'

deg_path <- v2_file('fig1_supp', 'GTE009_Sub1vs2_DEG.csv')
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

go <- enrich_go_dotplot(
  genes    = sub1_up,
  title    = 'Fig 1I — GO BP of Subclone_1 up genes',
  out_pdf  = v2_file('fig1_supp', 'Fig1I_GO_barplot.pdf'),
  csv_path = v2_file('fig1_supp', 'Fig1I_GO_subclone1.csv'),
  csv_writer = 'base',
  pcutoff = 0.05, qcutoff = 0.05
)
saveRDS(go, v2_file('fig1_supp', 'Fig1I_GO_subclone1.rds'))

cat('\n✓ Fig 1G/H/I (v2) done\n')
