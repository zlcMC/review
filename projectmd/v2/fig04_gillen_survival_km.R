#!/usr/bin/env Rscript
# fig04_gillen_survival_km.R (v2) — Gillen GSE125861 KM using Table S1 OS/PFS endpoints.
#
# v2 changes:
#   * Uses helpers::load_gse_expression + helpers::score_trajectory
#     to remove the duplicated GEO loader and trajectory scorer.
#   * Writes outputs to output/v2/fig4_survival/.

this_script <- sub('^--file=', '',
                   grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(this_script), 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(readxl)
  library(survival); library(survminer)
})

out       <- v2_dir('fig4_survival')
cache_dir <- output_path('external_cache')
pmc_dir   <- file.path(cache_dir, 'pmc_supp')
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

load_gse125861 <- function() {
  url  <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125861/matrix/GSE125861_series_matrix.txt.gz'
  dest <- file.path(cache_dir, 'GSE125861_series_matrix.txt.gz')
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    download.file(url, dest, mode = 'wb', quiet = FALSE)
  }
  dat <- load_gse_expression(
    matrix_path = dest,
    extra_columns = list(
      tumor_type = 'tumor type',
      tissue     = 'tissue'
    ),
    sample_filter_key    = 'tumor_type',
    sample_filter_values = c('primary', 'recurrent'),
    require_log2 = FALSE
  )
  dat$meta$sample_key <- sub('.*\\(([^)]+)\\).*', '\\1', dat$meta$title)
  dat
}

load_gillen_clinical <- function() {
  xlsx <- file.path(pmc_dir, 'Gillen_supplement_2.xlsx')
  if (!file.exists(xlsx)) {
    stop('Missing downloaded Gillen supplement: ', xlsx)
  }
  clinical <- read_excel(xlsx, sheet = 'Patient and sample details', skip = 1)
  pres_col <- 'presentation (1= presentation, 2= 1st recurrence)'
  clinical %>% transmute(
    sample_key   = paste0(as.integer(`patient ID`), '-', as.integer(.data[[pres_col]])),
    patient_id   = as.integer(`patient ID`),
    presentation = as.integer(.data[[pres_col]]),
    molecular_group_methyl = `methylomic molecular subgroup (using MolecularNeuropathology.org, https://www.molecularneuropathology.org/mnp)`,
    molecular_group_transcriptomic = `transcriptomic molecular subgroup (using non-negative matrix factorization (NMF), https://cloud.genepattern.org)`,
    who_grade        = `WHO grade`,
    age_years        = as.numeric(`age at Dx (years)`),
    sex              = sex,
    PFS_time_months  = as.numeric(`PFS (months)`),
    PFS_event        = as.integer(`recurrence event`),
    OS_time_months   = as.numeric(`OS (months)`),
    OS_event         = as.integer(`death event`)
  )
}

merge_survival_cohort <- function(scored, clinical) {
  cohort <- scored %>%
    left_join(clinical, by = 'sample_key') %>%
    filter(tumor_type == 'primary',
           !is.na(OS_time_months), !is.na(OS_event),
           !is.na(PFS_time_months), !is.na(PFS_event))
  median_score <- median(cohort$trajectory_score, na.rm = TRUE)
  cohort %>% mutate(
    trajectory_group = factor(
      ifelse(trajectory_score >= median_score, 'High', 'Low'),
      levels = c('Low', 'High')
    )
  )
}

endpoint_stats <- function(cohort, time_col, event_col, endpoint) {
  surv_formula <- as.formula(sprintf('Surv(%s, %s) ~ trajectory_group', time_col, event_col))
  diff      <- survdiff(surv_formula, data = cohort)
  p_logrank <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
  cox_group <- coxph(surv_formula, data = cohort)
  cox_cont  <- coxph(as.formula(sprintf('Surv(%s, %s) ~ trajectory_score',
                                        time_col, event_col)), data = cohort)

  data.frame(
    endpoint            = endpoint,
    n                   = nrow(cohort),
    events              = sum(cohort[[event_col]], na.rm = TRUE),
    median_score        = median(cohort$trajectory_score, na.rm = TRUE),
    low_n               = sum(cohort$trajectory_group == 'Low'),
    high_n              = sum(cohort$trajectory_group == 'High'),
    logrank_p           = p_logrank,
    group_HR_high_vs_low = unname(exp(coef(cox_group))[1]),
    group_HR_low95      = unname(exp(confint(cox_group))[1, 1]),
    group_HR_high95     = unname(exp(confint(cox_group))[1, 2]),
    group_p             = unname(summary(cox_group)$coefficients[1, 'Pr(>|z|)']),
    continuous_HR       = unname(exp(coef(cox_cont))[1]),
    continuous_p        = unname(summary(cox_cont)$coefficients[1, 'Pr(>|z|)']),
    stringsAsFactors    = FALSE
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
    fit, data = cohort, pval = FALSE, risk.table = TRUE, conf.int = FALSE,
    palette = c('#4DBBD5', '#E64B35'),
    legend.title = 'Trajectory score', legend.labs = c('Low', 'High'),
    title = paste('Fig. 4 survival:', endpoint, 'by trajectory score'),
    subtitle = subtitle, xlab = 'Months', ylab = paste(endpoint, 'probability'),
    risk.table.height = 0.24, ggtheme = theme_classic(base_size = 11)
  )
  pdf(file.path(out, paste0(filename_stub, '.pdf')), width = 7.2, height = 8.0)
  print(plot); dev.off()
  png(file.path(out, paste0(filename_stub, '.png')), width = 1500, height = 1650, res = 220)
  print(plot); dev.off()
  stats
}

dat      <- load_gse125861()
scored   <- score_trajectory(dat$expr, dat$meta)
clinical <- load_gillen_clinical()
cohort   <- merge_survival_cohort(scored, clinical)
if (nrow(cohort) < 10) {
  stop('Too few matched primary samples with OS/PFS endpoints: ', nrow(cohort))
}
write_csv(cohort, file.path(out, 'GSE125861_Gillen_survival_scored.csv'))

stats <- bind_rows(
  plot_endpoint(cohort, 'OS_time_months',  'OS_event',  'Overall survival',
                'Fig4_GSE125861_Gillen_trajectory_OS_KM'),
  plot_endpoint(cohort, 'PFS_time_months', 'PFS_event', 'Progression-free survival',
                'Fig4_GSE125861_Gillen_trajectory_PFS_KM')
)
stats$matched_undiff_genes <- attr(scored, 'matched_undiff')
stats$matched_diff_genes   <- attr(scored, 'matched_diff')
write_csv(stats, file.path(out, 'GSE125861_Gillen_survival_stats.csv'))

cat(sprintf('[v2] Matched %d primary GSE125861 samples with Gillen Table S1 OS/PFS.\n', nrow(cohort)))
cat(sprintf('[v2] Matched trajectory genes: undiff=%d, diff=%d.\n',
            attr(scored, 'matched_undiff'), attr(scored, 'matched_diff')))
cat('[v2] Wrote KM plots and stats to ', out, '\n', sep = '')
