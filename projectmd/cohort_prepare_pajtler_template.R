#!/usr/bin/env Rscript
# cohort_prepare_pajtler_template.R — 把 cBioPortal Pajtler 2015 EPN 数据
#   (epnt_mskcc_2024.tar.gz) 转成 fig04_survival_template.R 期望的格式
# 用法：
#   1) wget https://cbioportal-datahub.s3.amazonaws.com/epnt_mskcc_2024.tar.gz
#   2) mkdir -p projectfile/external_cohort/Pajtler2015
#   3) tar xzf epnt_mskcc_2024.tar.gz -C projectfile/external_cohort/Pajtler2015 --strip-components=1
#   4) conda run -n epn2_r Rscript projectmd/cohort_prepare_pajtler_template.R
#   5) conda run -n epn2_r Rscript projectmd/fig04_survival_template.R   # 自动出 KM+Cox

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', script_args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
    script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
    setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')

cdir <- raw_data_path('external_cohort/Pajtler2015')
stopifnot(dir.exists(cdir))

# 1) 表达矩阵
expr_files <- list.files(cdir, pattern='^data_mrna.*rsem.*\\.txt$|^data_mrna.*median\\.txt$', full.names=TRUE)
if (length(expr_files) == 0) stop('找不到 data_mrna_*.txt')
expr_f <- expr_files[1]
cat('读取', expr_f, '\n')
expr <- read.delim(expr_f, check.names=FALSE)
# 第一列 Hugo_Symbol，第二列 Entrez_Gene_Id（如有）；之后为样本
sym_col <- grep('Hugo|Symbol', colnames(expr), value=TRUE)[1]
keep_cols <- setdiff(colnames(expr), c('Hugo_Symbol','Entrez_Gene_Id'))
expr_mat <- aggregate(expr[, keep_cols], by=list(gene=expr[[sym_col]]),
                      FUN=mean, na.rm=TRUE)
rownames(expr_mat) <- expr_mat$gene; expr_mat$gene <- NULL
expr_mat <- as.matrix(expr_mat)
# log2 if not yet
if (max(expr_mat, na.rm=TRUE) > 50) expr_mat <- log2(expr_mat + 1)
write.csv(expr_mat, file.path(cdir, 'expr.csv'))
cat(sprintf('expr.csv: %d gene × %d sample\n', nrow(expr_mat), ncol(expr_mat)))

# 2) 临床
clin_p <- read.delim(file.path(cdir, 'data_clinical_patient.txt'),
                     comment.char='#', check.names=FALSE)
clin_s <- tryCatch(read.delim(file.path(cdir, 'data_clinical_sample.txt'),
                              comment.char='#', check.names=FALSE),
                   error=function(e) NULL)
# 提取生存
os_t <- grep('OS_MONTHS|OS_TIME', colnames(clin_p), value=TRUE)[1]
os_e <- grep('OS_STATUS|OS_EVENT', colnames(clin_p), value=TRUE)[1]
clin_out <- data.frame(
    sample = clin_p$PATIENT_ID,
    OS_time = as.numeric(clin_p[[os_t]]),
    OS_event = as.integer(grepl('1|DECEASED|DEAD', clin_p[[os_e]]))
)
# 如果 sample 列不在 expr 里，尝试用 sample-level 表
if (sum(clin_out$sample %in% colnames(expr_mat)) == 0 && !is.null(clin_s)) {
    cat('用 sample 表映射...\n')
    map <- setNames(clin_s$PATIENT_ID, clin_s$SAMPLE_ID)
    new_clin <- data.frame(
        sample = colnames(expr_mat),
        patient = map[colnames(expr_mat)]
    )
    new_clin <- merge(new_clin, clin_out, by.x='patient', by.y='sample')
    clin_out <- new_clin[, c('sample','OS_time','OS_event')]
}
clin_out <- clin_out[!is.na(clin_out$OS_time) & !is.na(clin_out$OS_event), ]
write.csv(clin_out, file.path(cdir, 'clin.csv'), row.names=FALSE)
cat(sprintf('clin.csv: %d 患者，事件 %d\n', nrow(clin_out), sum(clin_out$OS_event)))
cat('\n下一步: conda run -n epn2_r Rscript projectmd/fig04_survival_template.R\n')
