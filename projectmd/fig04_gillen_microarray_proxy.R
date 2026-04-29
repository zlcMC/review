#!/usr/bin/env Rscript
# Gillen et al. external microarray proxy for Fig. 4 recurrence-related trajectory score.
# GSE125861 contains a GPL570 expression matrix with primary/recurrent labels, but
# the GEO series matrix does not include OS/PFS fields. This script therefore
# generates a recurrence proxy, not a KM survival curve.

if (!file.exists('workspace_paths.R')) {
    stop('Run this script from the workspace root so workspace_paths.R is available.')
}
source('workspace_paths.R')

suppressPackageStartupMessages({
    library(GEOquery)
    library(limma)
    library(ggplot2)
    library(dplyr)
    library(readr)
})

out <- output_path('fig4_survival')
cache_dir <- output_path('external_cache')
dir.create(out, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

score_with_signatures <- function(expr, meta) {
    sigs <- readRDS(output_path('fig3/trajectory_signatures.rds'))
    undiff_sig <- if (!is.null(sigs$undiff_genes)) sigs$undiff_genes else sigs$undiff
    diff_sig <- if (!is.null(sigs$diff_genes)) sigs$diff_genes else sigs$diff
    undiff <- intersect(undiff_sig, rownames(expr))
    diff <- intersect(diff_sig, rownames(expr))
    if (length(undiff) < 5 || length(diff) < 5) {
        stop('Too few signature genes matched GSE125861 expression matrix.')
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

load_gse125861 <- function() {
    url <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125861/matrix/GSE125861_series_matrix.txt.gz'
    dest <- file.path(cache_dir, 'GSE125861_series_matrix.txt.gz')
    if (!file.exists(dest) || file.info(dest)$size == 0) {
        download.file(url, dest, mode = 'wb', quiet = FALSE)
    }
    gse <- getGEO(filename = dest, getGPL = FALSE)
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
        tumor_type = extract('tumor type'),
        tissue = extract('tissue'),
        stringsAsFactors = FALSE
    )
    keep <- !is.na(meta$tumor_type) & meta$tumor_type %in% c('primary', 'recurrent')
    expr <- expr[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
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

plot_gse125861 <- function(scored) {
    stats <- scored %>%
        group_by(tumor_type) %>%
        summarise(n = n(), mean = mean(trajectory_score), median = median(trajectory_score), .groups = 'drop')
    pval <- wilcox.test(trajectory_score ~ tumor_type, data = scored)$p.value
    readr::write_csv(stats, file.path(out, 'GSE125861_bulk_trajectory_score_group_stats.csv'))

    p <- ggplot(scored, aes(tumor_type, trajectory_score, fill = tumor_type)) +
        geom_violin(trim = FALSE, alpha = 0.75) +
        geom_boxplot(width = 0.16, fill = 'white', outlier.size = 0.7) +
        geom_jitter(width = 0.12, size = 0.8, alpha = 0.55) +
        scale_fill_manual(values = c(primary = '#4DBBD5', recurrent = '#E64B35')) +
        theme_classic(base_size = 11) +
        labs(
            title = 'Fig. 4 proxy: Gillen GSE125861 bulk trajectory score',
            subtitle = sprintf('Wilcoxon p = %.3g; matched undiff=%d, diff=%d genes',
                               pval, attr(scored, 'matched_undiff'), attr(scored, 'matched_diff')),
            x = NULL,
            y = 'mean(undiff genes) - mean(diff genes)'
        ) +
        theme(legend.position = 'none')
    ggsave(file.path(out, 'Fig4_GSE125861_bulk_trajectory_primary_vs_recurrent_proxy.pdf'), p,
           width = 5.2, height = 4.5)
    ggsave(file.path(out, 'Fig4_GSE125861_bulk_trajectory_primary_vs_recurrent_proxy.png'), p,
           width = 5.2, height = 4.5, dpi = 220)
}

plot_combined <- function(gse125861) {
    gse64415_path <- file.path(out, 'GSE64415_bulk_trajectory_scores.csv')
    if (!file.exists(gse64415_path)) return(invisible())
    gse64415 <- readr::read_csv(gse64415_path, show_col_types = FALSE) %>%
        transmute(dataset = 'GSE64415', condition = primary_relapse, trajectory_score)
    gillen <- gse125861 %>%
        transmute(dataset = 'GSE125861', condition = ifelse(tumor_type == 'recurrent', 'relapse/recurrent', 'primary'), trajectory_score)
    gse64415 <- gse64415 %>%
        mutate(condition = ifelse(condition == 'relapse', 'relapse/recurrent', condition))
    combined <- bind_rows(gse64415, gillen)
    readr::write_csv(combined, file.path(out, 'bulk_trajectory_proxy_two_cohorts.csv'))

    p <- ggplot(combined, aes(condition, trajectory_score, fill = condition)) +
        geom_boxplot(outlier.size = 0.45, alpha = 0.82) +
        geom_jitter(width = 0.12, size = 0.45, alpha = 0.35) +
        facet_wrap(~ dataset, scales = 'free_y') +
        scale_fill_manual(values = c(primary = '#4DBBD5', `relapse/recurrent` = '#E64B35')) +
        theme_bw(base_size = 9) +
        labs(title = 'Bulk trajectory score proxy in two external cohorts', x = NULL, y = 'trajectory score') +
        theme(legend.position = 'none', axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(file.path(out, 'Fig4_bulk_trajectory_two_cohorts_proxy.pdf'), p,
           width = 7.5, height = 4.2)
    ggsave(file.path(out, 'Fig4_bulk_trajectory_two_cohorts_proxy.png'), p,
           width = 7.5, height = 4.2, dpi = 220)
}

main <- function() {
    dat <- load_gse125861()
    scored <- score_with_signatures(dat$expr, dat$meta)
    readr::write_csv(scored, file.path(out, 'GSE125861_bulk_trajectory_scores.csv'))
    plot_gse125861(scored)
    plot_combined(scored)
    cat('Wrote GSE125861 recurrence proxy outputs to ', out, '\n', sep = '')
}

main()
