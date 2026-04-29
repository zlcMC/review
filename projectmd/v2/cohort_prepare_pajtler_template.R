#!/usr/bin/env Rscript
# v2 of cohort_prepare_pajtler_template.R — convert cBioPortal Pajtler 2015 dump
# into the format expected by fig04_survival_template.R.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

cdir <- raw_data_path('external_cohort/Pajtler2015')
stopifnot(dir.exists(cdir))

# 1) expression
expr_files <- list.files(cdir, full.names = TRUE,
  pattern = '^data_mrna.*rsem.*\\.txt$|^data_mrna.*median\\.txt$')
if (!length(expr_files)) stop('找不到 data_mrna_*.txt')
expr <- read.delim(expr_files[1], check.names = FALSE)
sym_col <- grep('Hugo|Symbol', colnames(expr), value = TRUE)[1]
keep    <- setdiff(colnames(expr), c('Hugo_Symbol', 'Entrez_Gene_Id'))
expr_mat <- aggregate(expr[, keep], by = list(gene = expr[[sym_col]]),
                      FUN = mean, na.rm = TRUE)
rownames(expr_mat) <- expr_mat$gene; expr_mat$gene <- NULL
expr_mat <- as.matrix(expr_mat)
if (max(expr_mat, na.rm = TRUE) > 50) expr_mat <- log2(expr_mat + 1)
write.csv(expr_mat, file.path(cdir, 'expr.csv'))
cat(sprintf('expr.csv: %d gene × %d sample\n', nrow(expr_mat), ncol(expr_mat)))

# 2) clinical
clin_p <- read.delim(file.path(cdir, 'data_clinical_patient.txt'),
                     comment.char = '#', check.names = FALSE)
clin_s <- tryCatch(read.delim(file.path(cdir, 'data_clinical_sample.txt'),
                              comment.char = '#', check.names = FALSE),
                   error = function(e) NULL)
os_t <- grep('OS_MONTHS|OS_TIME',  colnames(clin_p), value = TRUE)[1]
os_e <- grep('OS_STATUS|OS_EVENT', colnames(clin_p), value = TRUE)[1]
clin_out <- data.frame(
  sample   = clin_p$PATIENT_ID,
  OS_time  = as.numeric(clin_p[[os_t]]),
  OS_event = as.integer(grepl('1|DECEASED|DEAD', clin_p[[os_e]]))
)
if (sum(clin_out$sample %in% colnames(expr_mat)) == 0 && !is.null(clin_s)) {
  cat('mapping via sample table…\n')
  map <- setNames(clin_s$PATIENT_ID, clin_s$SAMPLE_ID)
  new_clin <- data.frame(sample = colnames(expr_mat),
                         patient = map[colnames(expr_mat)])
  new_clin <- merge(new_clin, clin_out, by.x = 'patient', by.y = 'sample')
  clin_out <- new_clin[, c('sample', 'OS_time', 'OS_event')]
}
clin_out <- clin_out[!is.na(clin_out$OS_time) & !is.na(clin_out$OS_event), ]
write.csv(clin_out, file.path(cdir, 'clin.csv'), row.names = FALSE)
cat(sprintf('clin.csv: %d patients, events=%d\n',
            nrow(clin_out), sum(clin_out$OS_event)))
cat('\nNext: Rscript projectmd/v2/fig04_survival_template.R\n')
