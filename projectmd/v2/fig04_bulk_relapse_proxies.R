#!/usr/bin/env Rscript
# fig04_bulk_relapse_proxies.R (v2) — bulk-data proxy panels for Fig. 4 / Fig. 5
# relapse analyses.
#
# v2 changes:
#   * Uses helpers::init_workspace + helpers::load_gse_expression
#     + helpers::score_trajectory to drop ~80 lines of duplicated boilerplate.
#   * Writes outputs to output/v2/fig4_survival/.

this_script <- sub('^--file=', '',
                   grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(this_script), 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
})

out <- v2_dir('fig4_survival')

# ---- 1) load GSE64415 + score trajectory ----
load_gse64415 <- function() {
  load_gse_expression(
    matrix_path = raw_data_path('bulk/GSE64415_series_matrix.txt.gz'),
    extra_columns = list(
      primary_relapse = 'primary/relapse',
      subgroup        = 'molecular subgroup',
      location        = 'brain location',
      age             = 'age at diagnosis',
      gender          = 'gender'
    ),
    sample_filter_key    = 'primary_relapse',
    sample_filter_values = c('primary', 'relapse')
  )
}

# ---- 2) plotting ----
plot_bulk_trajectory <- function(scored) {
  stats <- scored %>%
    group_by(primary_relapse) %>%
    summarise(n = n(), mean = mean(trajectory_score),
              median = median(trajectory_score), .groups = 'drop')
  pval <- wilcox.test(trajectory_score ~ primary_relapse, data = scored)$p.value
  write_csv(stats, file.path(out, 'GSE64415_bulk_trajectory_score_group_stats.csv'))

  p <- ggplot(scored, aes(primary_relapse, trajectory_score, fill = primary_relapse)) +
    geom_violin(trim = FALSE, alpha = 0.75) +
    geom_boxplot(width = 0.16, fill = 'white', outlier.size = 0.7) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.55) +
    scale_fill_manual(values = c(primary = '#4DBBD5', relapse = '#E64B35')) +
    theme_classic(base_size = 11) +
    labs(
      title = 'Fig. 4 proxy: bulk trajectory score in GSE64415',
      subtitle = sprintf('Wilcoxon p = %.3g; matched undiff=%d, diff=%d genes',
                         pval, attr(scored, 'matched_undiff'),
                         attr(scored, 'matched_diff')),
      x = NULL, y = 'mean(undiff genes) - mean(diff genes)'
    ) +
    theme(legend.position = 'none')
  for (ext in c('pdf', 'png')) {
    ggsave(file.path(out, paste0('Fig4_bulk_trajectory_primary_vs_relapse_proxy.', ext)),
           p, width = 5.2, height = 4.5,
           dpi = if (ext == 'png') 220 else 300)
  }

  p_sub <- ggplot(scored, aes(primary_relapse, trajectory_score, fill = primary_relapse)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    geom_jitter(width = 0.12, size = 0.5, alpha = 0.45) +
    facet_wrap(~ subgroup, scales = 'free_y') +
    scale_fill_manual(values = c(primary = '#4DBBD5', relapse = '#E64B35')) +
    theme_bw(base_size = 8) +
    labs(title = 'Bulk trajectory score by molecular subgroup',
         x = NULL, y = 'trajectory score') +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(out, 'Fig4_bulk_trajectory_by_subgroup_proxy.pdf'),
         p_sub, width = 10, height = 7)
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
  write_csv(lr, file.path(out, 'GSE64415_bulk_key_LR_expression.csv'))

  gene_stats <- bind_rows(lapply(present, function(gene) {
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
  })) %>% arrange(wilcox_p)
  write_csv(gene_stats, file.path(out, 'GSE64415_bulk_key_LR_stats.csv'))

  top_genes <- intersect(unique(c('MDK', 'NCL', 'PTN', 'PTPRZ1',
                                  head(gene_stats$gene, 8))),
                         present)
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
  for (ext in c('pdf', 'png')) {
    ggsave(file.path(out, paste0('Fig5DE_bulk_key_LR_primary_vs_relapse_proxy.', ext)),
           p_box, width = 9, height = 6.2,
           dpi = if (ext == 'png') 220 else 300)
  }

  mean_mat <- lr %>%
    group_by(primary_relapse) %>%
    summarise(across(all_of(present), \(x) mean(x, na.rm = TRUE)), .groups = 'drop') %>%
    pivot_longer(-primary_relapse, names_to = 'gene', values_to = 'mean_expression') %>%
    group_by(gene) %>%
    mutate(z = as.numeric(scale(mean_expression))) %>%
    ungroup()
  p_heat <- ggplot(mean_mat, aes(primary_relapse, reorder(gene, z), fill = z)) +
    geom_tile(color = 'white', linewidth = 0.2) +
    scale_fill_gradient2(low = '#2166AC', mid = 'white', high = '#B2182B',
                         midpoint = 0, name = 'row Z') +
    theme_minimal(base_size = 8) +
    labs(title = 'Mean expression of key CellChat ligand/receptor genes',
         x = NULL, y = NULL)
  ggsave(file.path(out, 'Fig5DE_bulk_key_LR_mean_heatmap_proxy.pdf'),
         p_heat, width = 4.8, height = 7)
}

# ---- 3) main ----
dat    <- load_gse64415()
scored <- score_trajectory(dat$expr, dat$meta)
write_csv(scored, file.path(out, 'GSE64415_bulk_trajectory_scores.csv'))
plot_bulk_trajectory(scored)
plot_lr_proxy(dat$expr, scored)
cat('[v2] Wrote bulk relapse proxy outputs to ', out, '\n', sep = '')
