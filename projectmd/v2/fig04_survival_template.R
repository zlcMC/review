#!/usr/bin/env Rscript
# fig04_survival_template.R (v2) — Figure 4C external KM survival template.
#
# v2 changes:
#   * Uses helpers::init_workspace + helpers::score_trajectory.
#   * Writes outputs to output/v2/fig4_survival/.
#   * Reads trajectory_signatures via dep_file (v2 first, then legacy).

this_script <- sub('^--file=', '',
                   grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(this_script), 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggplot2); library(dplyr)
})

SAMPLES <- c('GTE001', 'GTE002', 'GTE009', 'GTE012')
FIG6    <- v2_dir('fig4_survival')

# ---- 1) per-sample trajectory_score distribution ----
all_scores <- do.call(rbind, lapply(SAMPLES, function(s) {
  rds <- dep_file('fig3', paste0(s, '_scores.rds'))
  if (!file.exists(rds)) return(NULL)
  df <- readRDS(rds); df$sample <- s; df
}))
if (!is.null(all_scores) && nrow(all_scores) > 0) {
  p <- ggplot(all_scores, aes(x = sample, y = trajectory_score, fill = sample)) +
    geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, fill = 'white') +
    theme_classic() +
    ggtitle('Trajectory score distribution (4 samples)')
  ggsave(file.path(FIG6, 'trajectory_score_per_sample.pdf'), p, width = 6, height = 5)
  cat('Wrote:', file.path(FIG6, 'trajectory_score_per_sample.pdf'), '\n')
}

# ---- 2) optional external cohort KM (drop-in directory layout) ----
cohort_root <- raw_data_path('external_cohort')

run_cohort <- function(cdir) {
  cid    <- basename(cdir)
  expr_f <- file.path(cdir, 'expr.csv')
  clin_f <- file.path(cdir, 'clin.csv')
  if (!file.exists(expr_f) || !file.exists(clin_f)) {
    cat('  [skip]', cid, 'missing expr.csv or clin.csv\n'); return(invisible())
  }
  expr <- as.matrix(read.csv(expr_f, row.names = 1, check.names = FALSE))
  clin <- read.csv(clin_f)
  stopifnot(all(c('sample', 'OS_time', 'OS_event') %in% colnames(clin)))
  common <- intersect(clin$sample, colnames(expr))
  expr   <- expr[, common]; clin <- clin[match(common, clin$sample), ]

  # score_trajectory expects meta$sample as the lookup column.
  scored <- score_trajectory(expr, clin, sample_col = 'sample')
  cat(sprintf('  %s: %d patients; matched undiff=%d, diff=%d\n',
              cid, nrow(scored),
              attr(scored, 'matched_undiff'), attr(scored, 'matched_diff')))
  scored$group <- as.character(scored$trajectory_group)
  write.csv(scored, file.path(FIG6, sprintf('%s_scored.csv', cid)), row.names = FALSE)

  fit <- survfit(Surv(OS_time, OS_event) ~ group, data = scored)
  p_km <- ggsurvplot(fit, data = scored, pval = TRUE, risk.table = TRUE,
                     palette = c('firebrick', 'steelblue'),
                     title = sprintf('%s: trajectory High vs Low OS', cid))
  pdf(file.path(FIG6, sprintf('Fig4C_KM_%s.pdf', cid)), width = 7, height = 8)
  print(p_km); dev.off()

  cox <- coxph(Surv(OS_time, OS_event) ~ trajectory_score, data = scored)
  capture.output(summary(cox),
                 file = file.path(FIG6, sprintf('%s_cox.txt', cid)))
}

if (dir.exists(cohort_root)) {
  cohorts <- list.dirs(cohort_root, recursive = FALSE)
  if (length(cohorts) > 0) {
    cat('\nFound external cohorts:', paste(basename(cohorts), collapse = ', '), '\n')
    for (cd in cohorts) run_cohort(cd)
  } else {
    cat('\nprojectfile/external_cohort/ is empty. Layout expected:\n',
        '  projectfile/external_cohort/<cohort>/expr.csv  (gene x sample)\n',
        '  projectfile/external_cohort/<cohort>/clin.csv  (sample, OS_time, OS_event)\n',
        sep = '')
  }
} else {
  cat('\nNo projectfile/external_cohort/ found; strict KM survival needs an external cohort drop-in.\n')
}

cat('\n[v2] fig04_survival_template.R done →', FIG6, '\n')
