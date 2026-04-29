#!/usr/bin/env Rscript
# v2 of projectmd/fig02_infercnv_run_gte009.R
# 与原脚本逻辑完全一致；仅切换为 v2 路径 helper。
# 仍然是 HPC 级耗时任务（数小时 / 显存与多线程要求），本地可不重跑。
# 默认输出：output/v2/fig2_infercnv_strict/GTE009/
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(infercnv)
})
this_script <- sub('^--file=', '',
                   grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(this_script), "helpers", "helpers.R"))
init_workspace()

args <- commandArgs(trailingOnly = TRUE)
gene_order_file <- if (length(args) >= 1) args[[1]] else dep_file("reference", "gencode_hg38_gene_order.tsv")
out_dir         <- if (length(args) >= 2) args[[2]] else v2_dir("fig2_infercnv_strict", "GTE009")
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))

seurat_file <- raw_data_path("GTE009_seurat.rds")
stopifnot(file.exists(seurat_file), file.exists(gene_order_file))

input_dir <- file.path(out_dir, "input"); dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(seurat_file); DefaultAssay(obj) <- "RNA"
rna <- obj[["RNA"]]
counts <- if (inherits(rna, "Assay5")) {
  LayerData(rna, layer = "counts")
} else {
  GetAssayData(obj, assay = "RNA", slot = "counts")
}

meta <- obj@meta.data
stopifnot(all(c("CNV_Cluster", "cell_type") %in% colnames(meta)))

group <- rep(NA_character_, nrow(meta))
ref_ct <- c("Mic", "EC", "TC")
is_ref <- meta$cell_type %in% ref_ct
group[is_ref] <- paste0("Ref_", meta$cell_type[is_ref])
is_sub <- meta$CNV_Cluster %in% c("Subclone_1", "Subclone_2")
group[is_sub] <- as.character(meta$CNV_Cluster[is_sub])

keep <- !is.na(group)
counts <- counts[, keep, drop = FALSE]; group <- group[keep]

gene_order <- read.table(gene_order_file, sep = "\t", header = FALSE,
                         stringsAsFactors = FALSE,
                         col.names = c("gene", "chr", "start", "end"))
gene_order <- gene_order[!duplicated(gene_order$gene), ]
common <- intersect(rownames(counts), gene_order$gene)
stopifnot(length(common) >= 5000)
counts <- counts[common, , drop = FALSE]

ann_file <- file.path(input_dir, "GTE009_infercnv_annotations.tsv")
write.table(data.frame(cell = colnames(counts), group = group),
            ann_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(data.frame(
  metric = c("cells", "genes", "reference_groups", "tumor_groups", "threads"),
  value  = c(ncol(counts), nrow(counts),
             paste(sort(unique(group[grepl("^Ref_", group)])), collapse = ","),
             paste(sort(unique(group[!grepl("^Ref_", group)])), collapse = ","),
             threads)),
  file.path(input_dir, "GTE009_infercnv_input_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

message("Genes: ", nrow(counts), "  Cells: ", ncol(counts)); print(table(group))

icnv <- CreateInfercnvObject(
  raw_counts_matrix = counts, annotations_file = ann_file, delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = sort(unique(group[grepl("^Ref_", group)])))

icnv <- infercnv::run(icnv, cutoff = 0.1, out_dir = out_dir,
                      cluster_by_groups = TRUE, denoise = TRUE, HMM = FALSE,
                      num_threads = threads, plot_steps = FALSE,
                      no_prelim_plot = TRUE, no_plot = FALSE)
saveRDS(icnv, file.path(out_dir, "GTE009_infercnv_object.rds"))
message("inferCNV completed: ", out_dir)
