#!/usr/bin/env Rscript
# fig01_bulk_relapse_deg.R — Figure 1F: bulk EPN relapse vs primary DEG，与 sc DEG 对照
# 用 GSE64415 (Pajtler 2015, 209 EPN, 自带 primary/relapse 标签)
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', script_args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
    script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
    setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')
suppressPackageStartupMessages({
    library(GEOquery); library(limma); library(VennDiagram); library(grid)
    library(fgsea)
})
out <- output_path('fig1_bulk'); dir.create(out, showWarnings=FALSE, recursive=TRUE)

# ---- 1) load GSE64415 ----
f <- raw_data_path('bulk/GSE64415_series_matrix.txt.gz')
gse <- getGEO(filename=f, getGPL=FALSE)
expr <- exprs(gse); pd <- pData(gse)
cat(sprintf('expr: %d × %d\n', nrow(expr), ncol(expr)))

# ---- 1b) 解析 primary/relapse + subgroup ----
chars <- pd[, grep('characteristics', colnames(pd))]
extract <- function(key) {
    apply(chars, 1, function(r) {
        x <- grep(key, r, value=TRUE, ignore.case=TRUE)
        if (length(x)) sub('.*: *','',x[1]) else NA
    })
}
label   <- extract('primary/relapse')
subgrp  <- extract('molecular subgroup')
cat('primary/relapse:\n'); print(table(label, useNA='always'))
cat('subgroup:\n'); print(table(subgrp, useNA='always'))

keep <- !is.na(label) & label %in% c('primary','relapse')
expr <- expr[, keep]; group <- factor(label[keep], levels=c('primary','relapse'))
cat(sprintf('保留 %d 样本: primary=%d, relapse=%d\n',
            ncol(expr), sum(group=='primary'), sum(group=='relapse')))

if (max(expr, na.rm=TRUE) > 50) expr <- log2(expr + 1)
expr <- expr[rowSums(is.na(expr))==0, ]

# ---- 2) 探针 → gene symbol（GPL570 = hgu133plus2） ----
sample_id <- rownames(expr)[1]
cat('rowname 示例:', sample_id, '\n')
if (grepl('_at$', sample_id)) {
    suppressPackageStartupMessages(library(hgu133plus2.db))
    map <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(expr),
                                 columns='SYMBOL', keytype='PROBEID')
    map <- map[!is.na(map$SYMBOL) & !duplicated(map$PROBEID), ]
    expr <- expr[map$PROBEID, ]
    expr <- limma::avereps(expr, ID=map$SYMBOL)
    cat(sprintf('注释后 %d genes\n', nrow(expr)))
}

# ---- 3) limma DEG ----
design <- model.matrix(~ group)
fit <- eBayes(lmFit(expr, design))
deg <- topTable(fit, coef=2, number=Inf, adjust.method='BH')
deg$gene <- rownames(deg)
write.csv(deg, file.path(out,'bulk_DEG_relapse_vs_primary.csv'), row.names=FALSE)
cat('Top 10 bulk DEG (relapse vs primary):\n')
print(head(deg[, c('gene','logFC','P.Value','adj.P.Val')], 10))

# ---- 4) sc DEG 交集 ----
sc_path <- output_path('fig1_supp/GTE009_Sub1vs2_DEG.csv')
if (!file.exists(sc_path)) {
    cands <- list.files(output_path('fig1_supp'), pattern='DEG|marker', full.names=TRUE)
    if (length(cands)) sc_path <- cands[1]
}
if (file.exists(sc_path)) {
    sc <- read.csv(sc_path)
    cat(sprintf('sc DEG: %s (%d rows)\n', sc_path, nrow(sc)))
    fc_col <- grep('log2FC|logFC|avg_log', colnames(sc), value=TRUE, ignore.case=TRUE)[1]
    p_col  <- grep('p_val_adj|adj.P|padj', colnames(sc), value=TRUE, ignore.case=TRUE)[1]
    g_col  <- intersect(c('gene','symbol','X','Gene'), colnames(sc))[1]
    if (is.na(g_col)) g_col <- colnames(sc)[1]
    sc_up <- sc[[g_col]][sc[[fc_col]] >  0.25 & sc[[p_col]] < 0.05]
    sc_dn <- sc[[g_col]][sc[[fc_col]] < -0.25 & sc[[p_col]] < 0.05]
    bulk_up <- deg$gene[deg$logFC >  0.5 & deg$adj.P.Val < 0.1]
    bulk_dn <- deg$gene[deg$logFC < -0.5 & deg$adj.P.Val < 0.1]
    cat(sprintf('sc up=%d dn=%d  | bulk up=%d dn=%d\n',
                length(sc_up), length(sc_dn), length(bulk_up), length(bulk_dn)))
    ov_up <- intersect(sc_up, bulk_up); ov_dn <- intersect(sc_dn, bulk_dn)
    write.csv(data.frame(gene=ov_up), file.path(out,'sc_vs_bulk_overlap_up.csv'), row.names=FALSE)
    write.csv(data.frame(gene=ov_dn), file.path(out,'sc_vs_bulk_overlap_dn.csv'), row.names=FALSE)
    cat(sprintf('交集 up=%d dn=%d\n', length(ov_up), length(ov_dn)))

    pdf(file.path(out,'Fig1F_venn_up.pdf'), width=5, height=5)
    grid.draw(venn.diagram(list(sc_Sub1_up=sc_up, bulk_relapse_up=bulk_up),
                           filename=NULL, fill=c('#E64B35','#4DBBD5'), alpha=0.5))
    dev.off()
    pdf(file.path(out,'Fig1F_venn_dn.pdf'), width=5, height=5)
    grid.draw(venn.diagram(list(sc_Sub1_dn=sc_dn, bulk_relapse_dn=bulk_dn),
                           filename=NULL, fill=c('#3C5488','#00A087'), alpha=0.5))
    dev.off()

    ranks <- deg$logFC; names(ranks) <- deg$gene
    ranks <- sort(ranks[!is.na(ranks) & !duplicated(names(ranks))], decreasing=TRUE)
    pw <- list(sc_Sub1_up=sc_up, sc_Sub1_dn=sc_dn)
    pw <- pw[lengths(pw) >= 5]
    if (length(pw)) {
        res <- fgsea(pathways=pw, stats=ranks, minSize=5)
        res_out <- as.data.frame(res)
        res_out$leadingEdge <- sapply(res_out$leadingEdge, paste, collapse=';')
        write.csv(res_out, file.path(out,'GSEA_sc_in_bulk.csv'), row.names=FALSE)
        print(res[, c('pathway','pval','padj','NES','size')])
    }
} else {
    cat('未找到 sc DEG，跳过交集分析\n')
}
cat('\n✓ Fig 1F bulk DEG 完成 →', out, '\n')
