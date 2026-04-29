#!/usr/bin/env Rscript
# Bulk-data proxy panels for unfinished Fig. 4 / Fig. 5 relapse analyses.
# Run from workspace root:
#   conda run -n epn2_r Rscript projectmd/fig04_bulk_relapse_proxies.R

if (!file.exists('workspace_paths.R')) {
    stop('Run this script from the workspace root so workspace_paths.R is available.')
}
source('workspace_paths.R')

suppressPackageStartupMessages({
    library(GEOquery)
    library(limma)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(readr)
})

out <- output_path('fig4_survival')
dir.create(out, showWarnings = FALSE, recursive = TRUE)

load_gse64415 <- function() {
    gse <- getGEO(filename = raw_data_path('bulk/GSE64415_series_matrix.txt.gz'), getGPL = FALSE)
    expr <- exprs(gse)
    pd <- pData(gse)
    chars <- pd[, grep('characteristics', colnames(pd)), drop = FALSE]
    extract <- function(key) {
        apply(chars, 1, function(row) {
            hit <- grep(key, row, value = TRUE, ignore.case = TRUE)
            if (length(hit)) sub('.*: *', '', hit[1]) else NA_character_
        })
    }
    meta <- data.frame(
        sample = rownames(pd),
        title = pd$title,
        primary_relapse = extract('primary/relapse'),
        subgroup = extract('molecular subgroup'),
        location = extract('brain location'),
        age = suppressWarnings(as.numeric(extract('age at diagnosis'))),
        gender = extract('gender'),
        stringsAsFactors = FALSE
    )
    keep <- !is.na(meta$primary_relapse) & meta$primary_relapse %in% c('primary', 'relapse')
    expr <- expr[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
    if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)
    expr <- expr[rowSums(is.na(expr)) == 0, , drop = FALSE]

    if (grepl('_at$', rownames(expr)[1])) {
        suppressPackageStartupMessages(library(hgu133plus2.db))
        map <- suppressMessages(AnnotationDbi::select(
            hgu133plus2.db,
            keys = rownames(expr),
            columns = 'SYMBOL',
            keytype = 'PROBEID'
        ))
        map <- map[!is.na(map$SYMBOL) & !duplicated(map$PROBEID), ]
        expr <- expr[map$PROBEID, , drop = FALSE]
        expr <- limma::avereps(expr, ID = map$SYMBOL)
    }
    list(expr = expr, meta = meta)
}

score_bulk_trajectory <- function(expr, meta) {
    sigs <- readRDS(output_path('fig3/trajectory_signatures.rds'))
    undiff_sig <- if (!is.null(sigs$undiff_genes)) sigs$undiff_genes else sigs$undiff
    diff_sig <- if (!is.null(sigs$diff_genes)) sigs$diff_genes else sigs$diff
    undiff <- intersect(undiff_sig, rownames(expr))
    diff <- intersect(diff_sig, rownames(expr))
    if (length(undiff) < 5 || length(diff) < 5) {
        stop('Too few signature genes matched bulk expression matrix.')
    }
    score <- colMeans(expr[undiff, , drop = FALSE], na.rm = TRUE) -
        colMeans(expr[diff, , drop = FALSE], na.rm = TRUE)
    meta$trajectory_score <- as.numeric(score[meta$sample])
    meta$trajectory_group <- ifelse(
        meta$trajectory_score >= median(meta$trajectory_score, na.rm = TRUE),
        'High', 'Low'
    )
    attr(meta, 'matched_undiff') <- length(undiff)
    attr(meta, 'matched_diff') <- length(diff)
    meta
}

plot_bulk_trajectory <- function(scored) {
    stats <- scored %>%
        group_by(primary_relapse) %>%
        summarise(n = n(), mean = mean(trajectory_score), median = median(trajectory_score), .groups = 'drop')
    pval <- wilcox.test(trajectory_score ~ primary_relapse, data = scored)$p.value
    readr::write_csv(stats, file.path(out, 'GSE64415_bulk_trajectory_score_group_stats.csv'))

    p <- ggplot(scored, aes(primary_relapse, trajectory_score, fill = primary_relapse)) +
        geom_violin(trim = FALSE, alpha = 0.75) +
        geom_boxplot(width = 0.16, fill = 'white', outlier.size = 0.7) +
        geom_jitter(width = 0.12, size = 0.8, alpha = 0.55) +
        scale_fill_manual(values = c(primary = '#4DBBD5', relapse = '#E64B35')) +
        theme_classic(base_size = 11) +
        labs(
            title = 'Fig. 4 proxy: bulk trajectory score in GSE64415',
            subtitle = sprintf('Wilcoxon p = %.3g; matched undiff=%d, diff=%d genes',
                               pval, attr(scored, 'matched_undiff'), attr(scored, 'matched_diff')),
            x = NULL,
            y = 'mean(undiff genes) - mean(diff genes)'
        ) +
        theme(legend.position = 'none')
    ggsave(file.path(out, 'Fig4_bulk_trajectory_primary_vs_relapse_proxy.pdf'), p,
           width = 5.2, height = 4.5)
    ggsave(file.path(out, 'Fig4_bulk_trajectory_primary_vs_relapse_proxy.png'), p,
           width = 5.2, height = 4.5, dpi = 220)

    p_sub <- ggplot(scored, aes(primary_relapse, trajectory_score, fill = primary_relapse)) +
        geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
        geom_jitter(width = 0.12, size = 0.5, alpha = 0.45) +
        facet_wrap(~ subgroup, scales = 'free_y') +
        scale_fill_manual(values = c(primary = '#4DBBD5', relapse = '#E64B35')) +
        theme_bw(base_size = 8) +
        labs(
            title = 'Bulk trajectory score by molecular subgroup',
            x = NULL,
            y = 'trajectory score'
        ) +
        theme(legend.position = 'none', axis.text.x = element_text(angle = 35, hjust = 1))
    ggsave(file.path(out, 'Fig4_bulk_trajectory_by_subgroup_proxy.pdf'), p_sub,
           width = 10, height = 7)
}

plot_lr_proxy <- function(expr, scored) {
    lr_genes <- unique(c(
        'MDK', 'NCL', 'PTN', 'PTPRZ1', 'SDC2', 'SDC3', 'LRP1',
        'APP', 'CD74', 'CXCR4', 'TREM2', 'TYROBP', 'PPIA', 'BSG',
        'SPP1', 'CD44', 'HBEGF', 'EGFR', 'JAG1', 'NOTCH1', 'NOTCH3',
        'LGALS9', 'P4HB', 'COL1A1', 'COL1A2', 'COL4A1', 'COL6A1', 'COL6A2'
    ))
    present <- intersect(lr_genes, rownames(expr))
    lr <- t(expr[present, scored$sample, drop = FALSE]) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('sample') %>%
        left_join(scored, by = 'sample')
    readr::write_csv(lr, file.path(out, 'GSE64415_bulk_key_LR_expression.csv'))

    gene_stats <- lapply(present, function(gene) {
        form <- as.formula(sprintf('`%s` ~ primary_relapse', gene))
        pval <- wilcox.test(form, data = lr)$p.value
        data.frame(
            gene = gene,
            primary_mean = mean(lr[[gene]][lr$primary_relapse == 'primary'], na.rm = TRUE),
            relapse_mean = mean(lr[[gene]][lr$primary_relapse == 'relapse'], na.rm = TRUE),
            relapse_minus_primary = mean(lr[[gene]][lr$primary_relapse == 'relapse'], na.rm = TRUE) -
                mean(lr[[gene]][lr$primary_relapse == 'primary'], na.rm = TRUE),
            wilcox_p = pval
        )
    }) %>% bind_rows() %>% arrange(wilcox_p)
    readr::write_csv(gene_stats, file.path(out, 'GSE64415_bulk_key_LR_stats.csv'))

    top_genes <- unique(c('MDK', 'NCL', 'PTN', 'PTPRZ1', head(gene_stats$gene, 8)))
    top_genes <- intersect(top_genes, present)
    long <- lr %>%
        dplyr::select(sample, primary_relapse, all_of(top_genes)) %>%
        pivot_longer(all_of(top_genes), names_to = 'gene', values_to = 'expression')

    p_box <- ggplot(long, aes(primary_relapse, expression, fill = primary_relapse)) +
        geom_boxplot(outlier.size = 0.5, alpha = 0.82) +
        geom_jitter(width = 0.12, size = 0.45, alpha = 0.35) +
        facet_wrap(~ gene, scales = 'free_y', ncol = 4) +
        scale_fill_manual(values = c(primary = '#4DBBD5', relapse = '#E64B35')) +
        theme_bw(base_size = 8) +
        labs(title = 'Fig. 5D-E proxy: key CellChat ligand/receptor genes in GSE64415 bulk',
             x = NULL, y = 'log2 expression') +
        theme(legend.position = 'none', axis.text.x = element_text(angle = 35, hjust = 1))
    ggsave(file.path(out, 'Fig5DE_bulk_key_LR_primary_vs_relapse_proxy.pdf'), p_box,
           width = 9, height = 6.2)
    ggsave(file.path(out, 'Fig5DE_bulk_key_LR_primary_vs_relapse_proxy.png'), p_box,
           width = 9, height = 6.2, dpi = 220)

    mean_mat <- lr %>%
        group_by(primary_relapse) %>%
        summarise(across(all_of(present), \(x) mean(x, na.rm = TRUE)), .groups = 'drop') %>%
        pivot_longer(-primary_relapse, names_to = 'gene', values_to = 'mean_expression') %>%
        group_by(gene) %>%
        mutate(z = as.numeric(scale(mean_expression))) %>%
        ungroup()
    p_heat <- ggplot(mean_mat, aes(primary_relapse, reorder(gene, z), fill = z)) +
        geom_tile(color = 'white', linewidth = 0.2) +
        scale_fill_gradient2(low = '#2166AC', mid = 'white', high = '#B2182B', midpoint = 0,
                             name = 'row Z') +
        theme_minimal(base_size = 8) +
        labs(title = 'Mean expression of key CellChat ligand/receptor genes',
             x = NULL, y = NULL)
    ggsave(file.path(out, 'Fig5DE_bulk_key_LR_mean_heatmap_proxy.pdf'), p_heat,
           width = 4.8, height = 7)
}

main <- function() {
    dat <- load_gse64415()
    scored <- score_bulk_trajectory(dat$expr, dat$meta)
    readr::write_csv(scored, file.path(out, 'GSE64415_bulk_trajectory_scores.csv'))
    plot_bulk_trajectory(scored)
    plot_lr_proxy(dat$expr, scored)
    cat('Wrote bulk relapse proxy outputs to ', out, '\n', sep = '')
}

main()
