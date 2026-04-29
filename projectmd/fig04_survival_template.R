#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: R (epn2_r)
#     language: R
#     name: ir_epn2_r
# ---

# %% [markdown]
# fig04_survival_template.R — 对应论文 Figure 4C（生存分析）
# 论文做法：
#   - 用 trajectory_score 的均值把患者分为 high / low
#   - 在 ref (5,6) 发表的 bulk / 外部队列上做 Kaplan-Meier
#   - 本仓库目前 **没有** ref (5,6) 的外部数据，本脚本只给出：
#       1) 四个 EPN 样本细胞层面的 trajectory_score 分布 (补充结果)
#       2) 留给用户的 TODO：下载 ref 5 (Gillen 2020 Cancer Cell) 和
#          ref 6 (Gojo 2020 Cancer Cell) 补充表后，把患者 × score × survival
#          合并到 surv_df；脚本下段的 survminer::ggsurvplot 即可出图
# 运行：conda run -n epn2_r Rscript projectmd/fig04_survival_template.R

# %%
if (!file.exists('workspace_paths.R')) {
    stop('Run this script from the workspace root so workspace_paths.R is available.')
}
source('workspace_paths.R')

# %%
suppressPackageStartupMessages({
    library(survival); library(survminer); library(ggplot2); library(dplyr)
})

# %%
SAMPLES <- c('GTE001','GTE002','GTE009','GTE012')
FIG6 <- output_path('fig4_survival')
dir.create(FIG6, showWarnings = FALSE, recursive = TRUE)

# %%
# ---- 1) 合并四个样本的 trajectory_score 分布 ----
all_scores <- do.call(rbind, lapply(SAMPLES, function(s) {
    rds <- output_path(paste0('fig3/', s, '_scores.rds'))
    if (!file.exists(rds)) return(NULL)
    df <- readRDS(rds); df$sample <- s; df
}))
if (!is.null(all_scores) && nrow(all_scores) > 0) {
    p <- ggplot(all_scores, aes(x = sample, y = trajectory_score, fill = sample)) +
         geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, fill = 'white') +
         theme_classic() +
         ggtitle('Trajectory score distribution (4 samples)')
    ggsave(file.path(FIG6, 'trajectory_score_per_sample.pdf'),
           p, width = 6, height = 5)
    cat('已输出:', file.path(FIG6, 'trajectory_score_per_sample.pdf'), '\n')
}

# %%
# ---- 2) 外部生存队列：检测用户上传的文件 ----
# 用户准备：projectfile/external_cohort/{cohort_id}/
#   - expr.csv  (gene × sample，log2 TPM 或 normalized)
#   - clin.csv  (sample, OS_time, OS_event, [subtype])
# 脚本会自动用 trajectory_signatures.rds 的 diff/undiff 基因算每个患者的
# trajectory_score = mean(undiff) − mean(diff)，然后做 KM。

cohort_root <- raw_data_path('external_cohort')
sigs <- if (file.exists(output_path('fig3/trajectory_signatures.rds')))
            readRDS(output_path('fig3/trajectory_signatures.rds')) else NULL

run_cohort <- function(cdir) {
    cid <- basename(cdir)
    expr_f <- file.path(cdir, 'expr.csv')
    clin_f <- file.path(cdir, 'clin.csv')
    if (!file.exists(expr_f) || !file.exists(clin_f)) {
        cat('  [skip]', cid, '缺 expr.csv 或 clin.csv\n'); return(invisible())
    }
    expr <- as.matrix(read.csv(expr_f, row.names = 1, check.names = FALSE))
    clin <- read.csv(clin_f)
    stopifnot(all(c('sample','OS_time','OS_event') %in% colnames(clin)))
    common <- intersect(clin$sample, colnames(expr))
    expr <- expr[, common]; clin <- clin[match(common, clin$sample), ]
    diff_sig <- if (!is.null(sigs$diff_genes)) sigs$diff_genes else sigs$diff
    undiff_sig <- if (!is.null(sigs$undiff_genes)) sigs$undiff_genes else sigs$undiff
    diff_g  <- intersect(diff_sig,  rownames(expr))
    undiff_g<- intersect(undiff_sig,rownames(expr))
    cat(sprintf('  %s: %d 患者; diff %d, undiff %d 基因匹配\n',
                cid, length(common), length(diff_g), length(undiff_g)))
    score <- colMeans(expr[undiff_g, , drop=FALSE], na.rm=TRUE) -
             colMeans(expr[diff_g,   , drop=FALSE], na.rm=TRUE)
    clin$trajectory_score <- score
    clin$group <- ifelse(score >= median(score), 'High', 'Low')
    write.csv(clin, file.path(FIG6, sprintf('%s_scored.csv', cid)), row.names=FALSE)

    fit <- survfit(Surv(OS_time, OS_event) ~ group, data = clin)
    p_km <- ggsurvplot(fit, data = clin, pval = TRUE, risk.table = TRUE,
                       palette = c('firebrick','steelblue'),
                       title = sprintf('%s: trajectory High vs Low OS', cid))
    pdf(file.path(FIG6, sprintf('Fig4C_KM_%s.pdf', cid)), width = 7, height = 8)
    print(p_km); dev.off()

    cox <- coxph(Surv(OS_time, OS_event) ~ trajectory_score, data = clin)
    capture.output(summary(cox),
                   file = file.path(FIG6, sprintf('%s_cox.txt', cid)))
}

if (!is.null(sigs) && dir.exists(cohort_root)) {
    cohorts <- list.dirs(cohort_root, recursive = FALSE)
    if (length(cohorts) > 0) {
        cat('
检测到外部队列:', paste(basename(cohorts), collapse=', '), '\n')
        for (cd in cohorts) run_cohort(cd)
    } else {
        cat('
projectfile/external_cohort/ 为空。请按以下结构放数据：
  projectfile/external_cohort/Gillen2020/expr.csv  (gene × sample)
  projectfile/external_cohort/Gillen2020/clin.csv  (sample, OS_time, OS_event)
\n')
    }
} else {
    if (is.null(sigs)) {
        cat('
缺 output/fig3/trajectory_signatures.rds\n')
    } else {
        cat('
未发现 projectfile/external_cohort/ 目录；严格 KM 生存分析需先放入外部队列。\n')
    }
}

# %%
cat('
----- 外部生存分析模板（手动版） -----\n')

# %%
cat(paste(readLines(textConnection("
# surv_df <- read.csv(raw_data_path('external_EPN_cohort_with_scores.csv'))
# surv_df$group <- ifelse(surv_df$trajectory_score >= mean(surv_df$trajectory_score),
#                         'High', 'Low')
# fit <- survfit(Surv(OS_time, OS_event) ~ group, data = surv_df)
# p_km <- ggsurvplot(fit, data = surv_df, pval = TRUE, risk.table = TRUE,
#                    palette = c('firebrick','steelblue'),
#                    title = 'Trajectory-high vs -low OS')
# pdf(file.path(FIG6, 'Fig4C_KM_trajectory.pdf'), width = 7, height = 8)
# print(p_km); dev.off()
")), collapse = '\n'))
cat('\n\n✓ fig04_survival_template.R 占位完成（等外部数据）\n')
