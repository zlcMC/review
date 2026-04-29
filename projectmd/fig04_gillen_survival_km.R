#!/usr/bin/env Rscript
# Gillen et al. 2020 GSE125861 survival analysis using Table S1 OS/PFS endpoints.

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
    library(readxl)
    library(survival)
    library(survminer)
})

out <- output_path('fig4_survival')
cache_dir <- output_path('external_cache')
pmc_dir <- file.path(cache_dir, 'pmc_supp')
dir.create(out, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

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
    meta$sample_key <- sub('.*\\(([^)]+)\\).*', '\\1', meta$title)

    if (grepl('_at$', rownames(expr)[1])) {
        suppressPackageStartupMessages(library(hgu133plus2.db))
        probe_map <- suppressMessages(AnnotationDbi::select(
            hgu133plus2.db,
            keys = rownames(expr),
            columns = 'SYMBOL',
            keytype = 'PROBEID'
        ))
        probe_map <- probe_map[!is.na(probe_map$SYMBOL) & !duplicated(probe_map$PROBEID), ]
        expr <- expr[probe_map$PROBEID, , drop = FALSE]
        expr <- limma::avereps(expr, ID = probe_map$SYMBOL)
    }
    list(expr = expr, meta = meta)
}

score_with_signatures <- function(expr, meta) {
    sig_path <- output_path('fig3/trajectory_signatures.rds')
    if (!file.exists(sig_path)) {
        stop('Missing trajectory signatures: ', sig_path)
    }
    sigs <- readRDS(sig_path)
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
    attr(meta, 'matched_undiff') <- length(undiff)
    attr(meta, 'matched_diff') <- length(diff)
    meta
}

load_gillen_clinical <- function() {
    xlsx <- file.path(pmc_dir, 'Gillen_supplement_2.xlsx')
    if (!file.exists(xlsx)) {
        stop('Missing downloaded Gillen supplement: ', xlsx)
    }
    clinical <- readxl::read_excel(xlsx, sheet = 'Patient and sample details', skip = 1)
    presentation_col <- 'presentation (1= presentation, 2= 1st recurrence)'
    clinical %>%
        transmute(
            sample_key = paste0(as.integer(`patient ID`), '-', as.integer(.data[[presentation_col]])),
            patient_id = as.integer(`patient ID`),
            presentation = as.integer(.data[[presentation_col]]),
            molecular_group_methyl = `methylomic molecular subgroup (using MolecularNeuropathology.org, https://www.molecularneuropathology.org/mnp)`,
            molecular_group_transcriptomic = `transcriptomic molecular subgroup (using non-negative matrix factorization (NMF), https://cloud.genepattern.org)`,
            who_grade = `WHO grade`,
            age_years = as.numeric(`age at Dx (years)`),
            sex = sex,
            PFS_time_months = as.numeric(`PFS (months)`),
            PFS_event = as.integer(`recurrence event`),
            OS_time_months = as.numeric(`OS (months)`),
            OS_event = as.integer(`death event`)
        )
}

merge_survival_cohort <- function(scored, clinical) {
    cohort <- scored %>%
        left_join(clinical, by = 'sample_key') %>%
        filter(tumor_type == 'primary') %>%
        filter(!is.na(OS_time_months), !is.na(OS_event), !is.na(PFS_time_months), !is.na(PFS_event))

    median_score <- median(cohort$trajectory_score, na.rm = TRUE)
    cohort %>%
        mutate(
            trajectory_group = ifelse(trajectory_score >= median_score, 'High', 'Low'),
            trajectory_group = factor(trajectory_group, levels = c('Low', 'High'))
        )
}

endpoint_stats <- function(cohort, time_col, event_col, endpoint) {
    surv_formula <- as.formula(sprintf('Surv(%s, %s) ~ trajectory_group', time_col, event_col))
    diff <- survdiff(surv_formula, data = cohort)
    p_logrank <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
    cox_group <- coxph(surv_formula, data = cohort)
    cox_cont <- coxph(as.formula(sprintf('Surv(%s, %s) ~ trajectory_score', time_col, event_col)), data = cohort)

    data.frame(
        endpoint = endpoint,
        n = nrow(cohort),
        events = sum(cohort[[event_col]], na.rm = TRUE),
        median_score = median(cohort$trajectory_score, na.rm = TRUE),
        low_n = sum(cohort$trajectory_group == 'Low'),
        high_n = sum(cohort$trajectory_group == 'High'),
        logrank_p = p_logrank,
        group_HR_high_vs_low = unname(exp(coef(cox_group))[1]),
        group_HR_low95 = unname(exp(confint(cox_group))[1, 1]),
        group_HR_high95 = unname(exp(confint(cox_group))[1, 2]),
        group_p = unname(summary(cox_group)$coefficients[1, 'Pr(>|z|)']),
        continuous_HR = unname(exp(coef(cox_cont))[1]),
        continuous_p = unname(summary(cox_cont)$coefficients[1, 'Pr(>|z|)']),
        stringsAsFactors = FALSE
    )
}

plot_endpoint <- function(cohort, time_col, event_col, endpoint, filename_stub) {
    surv_formula <- as.formula(sprintf('Surv(%s, %s) ~ trajectory_group', time_col, event_col))
    fit <- survfit(surv_formula, data = cohort)
    fit$call$formula <- surv_formula
    stats <- endpoint_stats(cohort, time_col, event_col, endpoint)
    subtitle <- sprintf(
        'Gillen GSE125861 primary samples; n=%d, events=%d, log-rank p=%.3g',
        stats$n, stats$events, stats$logrank_p
    )
    plot <- ggsurvplot(
        fit,
        data = cohort,
        pval = FALSE,
        risk.table = TRUE,
        conf.int = FALSE,
        palette = c('#4DBBD5', '#E64B35'),
        legend.title = 'Trajectory score',
        legend.labs = c('Low', 'High'),
        title = paste('Fig. 4 survival:', endpoint, 'by trajectory score'),
        subtitle = subtitle,
        xlab = 'Months',
        ylab = paste(endpoint, 'probability'),
        risk.table.height = 0.24,
        ggtheme = theme_classic(base_size = 11)
    )

    pdf(file.path(out, paste0(filename_stub, '.pdf')), width = 7.2, height = 8.0)
    print(plot)
    dev.off()
    png(file.path(out, paste0(filename_stub, '.png')), width = 1500, height = 1650, res = 220)
    print(plot)
    dev.off()
    stats
}

main <- function() {
    dat <- load_gse125861()
    scored <- score_with_signatures(dat$expr, dat$meta)
    clinical <- load_gillen_clinical()
    cohort <- merge_survival_cohort(scored, clinical)

    if (nrow(cohort) < 10) {
        stop('Too few matched primary samples with OS/PFS endpoints: ', nrow(cohort))
    }

    readr::write_csv(cohort, file.path(out, 'GSE125861_Gillen_survival_scored.csv'))

    stats <- bind_rows(
        plot_endpoint(cohort, 'OS_time_months', 'OS_event', 'Overall survival',
                      'Fig4_GSE125861_Gillen_trajectory_OS_KM'),
        plot_endpoint(cohort, 'PFS_time_months', 'PFS_event', 'Progression-free survival',
                      'Fig4_GSE125861_Gillen_trajectory_PFS_KM')
    )
    stats$matched_undiff_genes <- attr(scored, 'matched_undiff')
    stats$matched_diff_genes <- attr(scored, 'matched_diff')
    readr::write_csv(stats, file.path(out, 'GSE125861_Gillen_survival_stats.csv'))

    cat(sprintf('Matched %d primary GSE125861 samples with Gillen Table S1 OS/PFS endpoints.\n', nrow(cohort)))
    cat(sprintf('Matched trajectory genes: undiff=%d, diff=%d.\n', attr(scored, 'matched_undiff'), attr(scored, 'matched_diff')))
    cat('Wrote KM plots and stats to ', out, '\n', sep = '')
}

main()