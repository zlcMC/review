#!/usr/bin/env Rscript
# fig04_gillen_microarray_proxy.R (v2) — Gillen GSE125861 recurrence proxy.
#
# v2 changes:
#   * Uses helpers::load_gse_expression + helpers::score_trajectory
#     instead of duplicated load_gse125861 / score_with_signatures blocks.
#   * Writes outputs to output/v2/fig4_survival/.

this_script <- sub('^--file=', '',
                   grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(this_script), 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
})

out       <- v2_dir('fig4_survival')
cache_dir <- output_path('external_cache')
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

load_gse125861 <- function() {
  url  <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125861/matrix/GSE125861_series_matrix.txt.gz'
  dest <- file.path(cache_dir, 'GSE125861_series_matrix.txt.gz')
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    download.file(url, dest, mode = 'wb', quiet = FALSE)
  }
  load_gse_expression(
    matrix_path = dest,
    extra_columns = list(
      tumor_type = 'tumor type',
      tissue     = 'tissue'
    ),
    sample_filter_key    = 'tumor_type',
    sample_filter_values = c('primary', 'recurrent'),
    require_log2 = FALSE
  )
}

plot_gse125861 <- function(scored) {
  stats <- scored %>%
    group_by(tumor_type) %>%
    summarise(n = n(), mean = mean(trajectory_score),
              median = median(trajectory_score), .groups = 'drop')
  pval <- wilcox.test(trajectory_score ~ tumor_type, data = scored)$p.value
  write_csv(stats, file.path(out, 'GSE125861_bulk_trajectory_score_group_stats.csv'))

  p <- ggplot(scored, aes(tumor_type, trajectory_score, fill = tumor_type)) +
    geom_violin(trim = FALSE, alpha = 0.75) +
    geom_boxplot(width = 0.16, fill = 'white', outlier.size = 0.7) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.55) +
    scale_fill_manual(values = c(primary = '#4DBBD5', recurrent = '#E64B35')) +
    theme_classic(base_size = 11) +
    labs(
      title = 'Fig. 4 proxy: Gillen GSE125861 bulk trajectory score',
      subtitle = sprintf('Wilcoxon p = %.3g; matched undiff=%d, diff=%d genes',
                         pval, attr(scored, 'matched_undiff'),
                         attr(scored, 'matched_diff')),
      x = NULL, y = 'mean(undiff genes) - mean(diff genes)'
    ) +
    theme(legend.position = 'none')
  for (ext in c('pdf', 'png')) {
    ggsave(file.path(out, paste0('Fig4_GSE125861_bulk_trajectory_primary_vs_recurrent_proxy.', ext)),
           p, width = 5.2, height = 4.5,
           dpi = if (ext == 'png') 220 else 300)
  }
}

plot_combined <- function(gillen_scored) {
  gse64415_path <- dep_file('fig4_survival', 'GSE64415_bulk_trajectory_scores.csv')
  if (!file.exists(gse64415_path)) return(invisible())
  gse64415 <- read_csv(gse64415_path, show_col_types = FALSE) %>%
    transmute(dataset = 'GSE64415',
              condition = ifelse(primary_relapse == 'relapse',
                                 'relapse/recurrent', primary_relapse),
              trajectory_score)
  gillen <- gillen_scored %>%
    transmute(dataset = 'GSE125861',
              condition = ifelse(tumor_type == 'recurrent',
                                 'relapse/recurrent', 'primary'),
              trajectory_score)
  combined <- bind_rows(gse64415, gillen)
  write_csv(combined, file.path(out, 'bulk_trajectory_proxy_two_cohorts.csv'))

  p <- ggplot(combined, aes(condition, trajectory_score, fill = condition)) +
    geom_boxplot(outlier.size = 0.45, alpha = 0.82) +
    geom_jitter(width = 0.12, size = 0.45, alpha = 0.35) +
    facet_wrap(~ dataset, scales = 'free_y') +
    scale_fill_manual(values = c(primary = '#4DBBD5', `relapse/recurrent` = '#E64B35')) +
    theme_bw(base_size = 9) +
    labs(title = 'Bulk trajectory score proxy in two external cohorts',
         x = NULL, y = 'trajectory score') +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 25, hjust = 1))
  for (ext in c('pdf', 'png')) {
    ggsave(file.path(out, paste0('Fig4_bulk_trajectory_two_cohorts_proxy.', ext)),
           p, width = 7.5, height = 4.2,
           dpi = if (ext == 'png') 220 else 300)
  }
}

dat    <- load_gse125861()
scored <- score_trajectory(dat$expr, dat$meta)
write_csv(scored, file.path(out, 'GSE125861_bulk_trajectory_scores.csv'))
plot_gse125861(scored)
plot_combined(scored)
cat('[v2] Wrote GSE125861 recurrence proxy outputs to ', out, '\n', sep = '')
