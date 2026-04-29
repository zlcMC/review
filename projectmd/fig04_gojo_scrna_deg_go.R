#!/usr/bin/env Rscript
# GO enrichment for Gojo GSE141460 recurrent-high malignant pseudobulk genes.

if (!file.exists('workspace_paths.R')) {
    stop('Run this script from the workspace root so workspace_paths.R is available.')
}
source('workspace_paths.R')

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(readr)
    library(dplyr)
})

out <- output_path('fig4_gojo_scrna')
deg_path <- file.path(out, 'GSE141460_malignant_pseudobulk_DEG_PFA_only.csv')
if (!file.exists(deg_path)) {
    stop('Missing DEG file. Run projectmd/fig04_gojo_scrna_validation.py first: ', deg_path)
}

deg <- readr::read_csv(deg_path, show_col_types = FALSE)
deg <- deg %>% filter(!is.na(ttest_p), recurrent_minus_primary > 0)

candidate <- deg %>% filter(ttest_p < 0.05, recurrent_minus_primary > 0.25)
if (nrow(candidate) < 10) {
    candidate <- deg %>% arrange(ttest_p, desc(recurrent_minus_primary)) %>% head(150)
}
genes <- unique(candidate$gene)
readr::write_lines(genes, file.path(out, 'GSE141460_recurrent_high_GO_input_genes.txt'))

ego <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.2,
    readable = TRUE
)

ego_df <- as.data.frame(ego)
readr::write_csv(ego_df, file.path(out, 'GSE141460_PFA_recurrent_high_GO_BP.csv'))

if (nrow(ego_df) > 0) {
    p <- dotplot(ego, showCategory = min(15, nrow(ego_df))) +
        ggtitle('Gojo GSE141460 PF-A recurrent-high pseudobulk GO') +
        theme_bw(base_size = 9)
    ggsave(file.path(out, 'Fig4D_GSE141460_PFA_recurrent_high_GO_BP_scrna.pdf'), p, width = 7.4, height = 5.4)
    ggsave(file.path(out, 'Fig4D_GSE141460_PFA_recurrent_high_GO_BP_scrna.png'), p, width = 7.4, height = 5.4, dpi = 220)
} else {
    cat('No enriched GO terms passed the configured cutoffs.\n')
}

cat(sprintf('GO input genes: %d; enriched terms: %d\n', length(genes), nrow(ego_df)))