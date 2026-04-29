#!/usr/bin/env Rscript
# v2 add-on: Pajtler 2015 GSE64415 microarray cohort
#
# Goal: provide subtype-balanced PF_EPN_A primary vs relapse trajectory-score
#       comparison (Wu cohort all PF-A primary; Gillen relapse all ST-RELA;
#       Gojo PF-A relapse n=8 vs primary n=4). Pajtler PF-A: 56 primary + 8 relapse.
#
# Inputs (must already be downloaded):
#   output/external_cache/GSE64415_series_matrix.txt.gz
#   output/fig3/undiff_genes.txt
#   output/fig3/diff_genes.txt
#
# Outputs (under output/v2/fig4_external_pajtler/):
#   GSE64415_sample_metadata.csv
#   GSE64415_trajectory_scores.csv
#   GSE64415_subtype_trajectory_summary.csv
#   GSE64415_PFA_primary_vs_relapse_stats.csv
#   Fig4_pajtler_PFA_primary_vs_relapse.pdf/.png
#   Fig4_pajtler_subtype_trajectory_overview.pdf/.png

this_script <- sub('^--file=', '',
                   grep('^--file=', commandArgs(trailingOnly=FALSE), value=TRUE)[1])
source(file.path(dirname(this_script), 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(ggplot2)
  library(dplyr)
})

OUT_DIR <- v2_dir('fig4_external_pajtler')
SM_PATH <- file.path(output_dir, 'external_cache', 'GSE64415_series_matrix.txt.gz')
stopifnot(file.exists(SM_PATH))

cat('[1/6] Loading GSE64415 series matrix...\n')
gse <- getGEO(filename = SM_PATH, GSEMatrix = TRUE, getGPL = FALSE, AnnotGPL = FALSE)
expr <- exprs(gse)
pdat <- pData(gse)
cat(sprintf('  Probes x samples: %d x %d\n', nrow(expr), ncol(expr)))

# --- 2) Build clean sample metadata
cat('[2/6] Parsing sample metadata...\n')
parse_char <- function(field_value, key) {
  v <- as.character(field_value)
  out <- ifelse(grepl(paste0('^', key, ':\\s*'), v),
                sub(paste0('^', key, ':\\s*'), '', v),
                NA_character_)
  out
}
char_cols <- grep('^characteristics_ch1', colnames(pdat), value = TRUE)
char_df <- pdat[, char_cols, drop = FALSE]

extract_field <- function(key) {
  apply(char_df, 1, function(row) {
    hits <- parse_char(row, key)
    val <- hits[!is.na(hits)]
    if (length(val) == 0) NA_character_ else val[1]
  })
}

meta <- data.frame(
  gsm        = pdat$geo_accession,
  title      = pdat$title,
  subgroup   = extract_field('molecular subgroup'),
  condition  = extract_field('primary/relapse'),
  location   = extract_field('brain location'),
  age_years  = suppressWarnings(as.numeric(extract_field('age at diagnosis \\(years\\)'))),
  sex        = extract_field('gender'),
  stringsAsFactors = FALSE
)
write.csv(meta, file.path(OUT_DIR, 'GSE64415_sample_metadata.csv'), row.names = FALSE)
cat(sprintf('  Wrote sample metadata: %d samples\n', nrow(meta)))
cat(sprintf('  PF_EPN_A primary=%d relapse=%d\n',
            sum(meta$subgroup=='PF_EPN_A' & meta$condition=='primary', na.rm=TRUE),
            sum(meta$subgroup=='PF_EPN_A' & meta$condition=='relapse', na.rm=TRUE)))

# --- 3) Probe -> gene symbol mapping (best-probe by max IQR per gene)
cat('[3/6] Mapping probes to gene symbols...\n')
probe2sym <- AnnotationDbi::select(hgu133plus2.db,
                                   keys = rownames(expr),
                                   columns = c('SYMBOL'),
                                   keytype = 'PROBEID')
probe2sym <- probe2sym[!is.na(probe2sym$SYMBOL) & probe2sym$SYMBOL != '', ]
# unique probe (some probes map to multiple symbols -> pick first)
probe2sym <- probe2sym[!duplicated(probe2sym$PROBEID), ]
cat(sprintf('  %d probes have a SYMBOL annotation\n', nrow(probe2sym)))

expr_sub <- expr[probe2sym$PROBEID, , drop = FALSE]
sym_vec  <- probe2sym$SYMBOL
# pick probe with largest IQR per gene
iqr_vals <- apply(expr_sub, 1, IQR, na.rm = TRUE)
ord <- order(sym_vec, -iqr_vals)
expr_sub <- expr_sub[ord, , drop = FALSE]
sym_vec  <- sym_vec[ord]
keep <- !duplicated(sym_vec)
expr_gene <- expr_sub[keep, , drop = FALSE]
rownames(expr_gene) <- sym_vec[keep]
cat(sprintf('  Gene-level matrix: %d genes x %d samples\n',
            nrow(expr_gene), ncol(expr_gene)))

# --- 4) Trajectory score per sample
cat('[4/6] Computing trajectory score (undiff - diff)...\n')
undiff_genes <- readLines(dep_file('fig3', 'undiff_genes.txt'))
diff_genes   <- readLines(dep_file('fig3', 'diff_genes.txt'))
undiff_in <- intersect(undiff_genes, rownames(expr_gene))
diff_in   <- intersect(diff_genes,   rownames(expr_gene))
cat(sprintf('  undiff matched: %d / %d ; diff matched: %d / %d\n',
            length(undiff_in), length(undiff_genes),
            length(diff_in), length(diff_genes)))

score_undiff <- colMeans(expr_gene[undiff_in, , drop = FALSE], na.rm = TRUE)
score_diff   <- colMeans(expr_gene[diff_in,   , drop = FALSE], na.rm = TRUE)
trajectory_score <- score_undiff - score_diff

scores <- data.frame(
  gsm = colnames(expr_gene),
  undiff_score = score_undiff,
  diff_score   = score_diff,
  trajectory_score = trajectory_score,
  stringsAsFactors = FALSE
)
scores <- merge(meta, scores, by = 'gsm', sort = FALSE)
write.csv(scores, file.path(OUT_DIR, 'GSE64415_trajectory_scores.csv'), row.names = FALSE)

# --- 5) Subtype-level summary
cat('[5/6] Subtype-level summary...\n')
summary_df <- scores %>%
  group_by(subgroup) %>%
  summarise(n = n(),
            mean_traj = mean(trajectory_score),
            sd_traj   = sd(trajectory_score),
            median_traj = median(trajectory_score),
            .groups = 'drop')
write.csv(summary_df, file.path(OUT_DIR, 'GSE64415_subtype_trajectory_summary.csv'),
          row.names = FALSE)
print(summary_df)

# --- 6) PF_EPN_A primary vs relapse: Wilcoxon + multivariate linear model
cat('[6/6] PF-A primary vs relapse stats and plots...\n')
pfa <- scores[scores$subgroup == 'PF_EPN_A' & scores$condition %in% c('primary', 'relapse'), ]
pfa$condition <- factor(pfa$condition, levels = c('primary', 'relapse'))
n_p <- sum(pfa$condition == 'primary')
n_r <- sum(pfa$condition == 'relapse')
cat(sprintf('  PF_EPN_A: primary=%d, relapse=%d\n', n_p, n_r))

w <- wilcox.test(trajectory_score ~ condition, data = pfa, exact = FALSE)
t_obj <- t.test(trajectory_score ~ condition, data = pfa)
# linear model adjusted for age + sex
mdl_uni <- lm(trajectory_score ~ condition, data = pfa)
mdl_adj <- lm(trajectory_score ~ condition + age_years + sex, data = pfa)

stats_out <- data.frame(
  test = c('wilcoxon', 't_test', 'lm_univariate', 'lm_adj_age_sex'),
  estimate = c(NA,
               diff(rev(t_obj$estimate)),
               coef(mdl_uni)['conditionrelapse'],
               coef(mdl_adj)['conditionrelapse']),
  p_value = c(w$p.value, t_obj$p.value,
              summary(mdl_uni)$coefficients['conditionrelapse', 'Pr(>|t|)'],
              summary(mdl_adj)$coefficients['conditionrelapse', 'Pr(>|t|)']),
  n_primary = n_p, n_relapse = n_r,
  stringsAsFactors = FALSE
)
write.csv(stats_out, file.path(OUT_DIR, 'GSE64415_PFA_primary_vs_relapse_stats.csv'),
          row.names = FALSE)
print(stats_out)

# Plot 1: PF-A primary vs relapse
pdf_png <- function(p, base, w = 4.5, h = 4.5) {
  ggsave(paste0(base, '.pdf'), p, width = w, height = h)
  ggsave(paste0(base, '.png'), p, width = w, height = h, dpi = 300)
}

p1 <- ggplot(pfa, aes(x = condition, y = trajectory_score, fill = condition)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.85) +
  scale_fill_manual(values = c(primary = '#4C72B0', relapse = '#C44E52')) +
  labs(x = NULL, y = 'Trajectory score (undiff - diff)',
       title = sprintf('Pajtler GSE64415 PF_EPN_A\nWilcoxon p=%.3g (n=%d vs %d)',
                       w$p.value, n_p, n_r)) +
  theme_classic(base_size = 12) + theme(legend.position = 'none')
pdf_png(p1, file.path(OUT_DIR, 'Fig4_pajtler_PFA_primary_vs_relapse'), 4.6, 4.5)

# Plot 2: subtype overview (all 8 subtypes)
sub_order <- summary_df$subgroup[order(-summary_df$mean_traj)]
scores$subgroup <- factor(scores$subgroup, levels = sub_order)
p2 <- ggplot(scores, aes(x = subgroup, y = trajectory_score, fill = subgroup)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.0, alpha = 0.6) +
  labs(x = NULL, y = 'Trajectory score (undiff - diff)',
       title = 'Pajtler GSE64415 trajectory score by molecular subgroup (n=209)') +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = 'none')
pdf_png(p2, file.path(OUT_DIR, 'Fig4_pajtler_subtype_trajectory_overview'), 7.0, 4.6)

cat('\nDone. Outputs in: ', OUT_DIR, '\n', sep = '')
