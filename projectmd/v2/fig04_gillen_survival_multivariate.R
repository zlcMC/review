#!/usr/bin/env Rscript
# task 4 — Fig 4A/C multivariate Cox sensitivity analysis.
# 在 Gillen GSE125861 + Table S1 cohort 上加 multivariate Cox：
#   trajectory + age + sex + subtype (methyl) + WHO grade
# 输出 output/v2/fig4_survival/multivariate_cox_summary.csv

this_script <- sub('^--file=', '',
                   grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(this_script), 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(readxl); library(survival)
})

out <- v2_dir('fig4_survival')

# 复用 v1 的 KM script 已生成的 cohort 表
cohort_csv <- dep_file('fig4_survival', 'GSE125861_Gillen_survival_scored.csv')
if (!file.exists(cohort_csv)) {
  stop('Missing ', cohort_csv, '；请先跑 fig04_gillen_survival_km.R')
}
cohort <- read_csv(cohort_csv, show_col_types = FALSE)
cat('Loaded cohort: n =', nrow(cohort), '\n')
cat('Columns:', paste(colnames(cohort), collapse=', '), '\n\n')

# subtype 简化：把 methyl group 折叠成 PF-A / ST-RELA / Other
cohort <- cohort %>% mutate(
  subtype_group = case_when(
    grepl('PFA',  molecular_group_methyl, ignore.case = TRUE) ~ 'PF-A',
    grepl('PFB',  molecular_group_methyl, ignore.case = TRUE) ~ 'PF-B',
    grepl('RELA', molecular_group_methyl, ignore.case = TRUE) ~ 'ST-RELA',
    grepl('YAP',  molecular_group_methyl, ignore.case = TRUE) ~ 'ST-YAP',
    TRUE ~ 'Other'
  ),
  sex_male = ifelse(tolower(as.character(sex)) %in% c('m','male'), 1, 0)
)
cat('Subtype distribution:\n'); print(table(cohort$subtype_group, useNA = 'ifany'))

run_models <- function(time_col, event_col, endpoint) {
  d <- cohort %>% filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]),
                         !is.na(trajectory_score), !is.na(age_years), !is.na(sex_male),
                         !is.na(subtype_group))
  fmls <- list(
    univariate_continuous   = sprintf('Surv(%s, %s) ~ trajectory_score', time_col, event_col),
    plus_age_sex            = sprintf('Surv(%s, %s) ~ trajectory_score + age_years + sex_male', time_col, event_col),
    plus_age_sex_subtype    = sprintf('Surv(%s, %s) ~ trajectory_score + age_years + sex_male + subtype_group', time_col, event_col),
    plus_age_sex_subtype_grade = sprintf('Surv(%s, %s) ~ trajectory_score + age_years + sex_male + subtype_group + who_grade', time_col, event_col)
  )
  rows <- list()
  for (nm in names(fmls)) {
    m <- tryCatch(coxph(as.formula(fmls[[nm]]), data = d), error = function(e) NULL)
    if (is.null(m)) next
    s <- summary(m)
    if (!'trajectory_score' %in% rownames(s$coefficients)) next
    rows[[length(rows)+1]] <- data.frame(
      endpoint   = endpoint,
      model      = nm,
      n          = m$n,
      n_events   = m$nevent,
      trajectory_HR = unname(s$coefficients['trajectory_score', 'exp(coef)']),
      trajectory_HR_low95 = unname(s$conf.int['trajectory_score', 'lower .95']),
      trajectory_HR_high95 = unname(s$conf.int['trajectory_score', 'upper .95']),
      trajectory_p = unname(s$coefficients['trajectory_score', 'Pr(>|z|)']),
      formula     = fmls[[nm]]
    )
  }
  do.call(rbind, rows)
}

res_os  <- run_models('OS_time_months',  'OS_event',  'OS')
res_pfs <- run_models('PFS_time_months', 'PFS_event', 'PFS')
res <- rbind(res_os, res_pfs)
write.csv(res, file.path(out, 'multivariate_cox_summary.csv'), row.names = FALSE)
print(res)
cat(sprintf('\n→ %s/multivariate_cox_summary.csv\n', out))
