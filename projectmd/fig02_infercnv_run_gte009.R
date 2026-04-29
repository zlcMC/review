#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(infercnv)
})

source("workspace_paths.R")

args <- commandArgs(trailingOnly = TRUE)
gene_order_file <- if (length(args) >= 1) args[[1]] else output_path("reference", "gencode_hg38_gene_order.tsv")
out_dir <- if (length(args) >= 2) args[[2]] else output_path("fig2_infercnv_strict", "GTE009")
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))

seurat_file <- raw_data_path("GTE009_seurat.rds")
if (!file.exists(seurat_file)) {
  stop("Missing Seurat input: ", seurat_file)
}
if (!file.exists(gene_order_file)) {
  stop("Missing gene-order file: ", gene_order_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
input_dir <- file.path(out_dir, "input")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(seurat_file)
DefaultAssay(obj) <- "RNA"
rna <- obj[["RNA"]]

counts <- if (inherits(rna, "Assay5")) {
  LayerData(rna, layer = "counts")
} else {
  GetAssayData(obj, assay = "RNA", slot = "counts")
}

meta <- obj@meta.data
required_columns <- c("CNV_Cluster", "cell_type")
missing_columns <- setdiff(required_columns, colnames(meta))
if (length(missing_columns) > 0) {
  stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "))
}

group <- rep(NA_character_, nrow(meta))
reference_celltypes <- c("Mic", "EC", "TC")
is_reference <- meta$cell_type %in% reference_celltypes
group[is_reference] <- paste0("Ref_", meta$cell_type[is_reference])

is_subclone <- meta$CNV_Cluster %in% c("Subclone_1", "Subclone_2")
group[is_subclone] <- as.character(meta$CNV_Cluster[is_subclone])

keep_cells <- !is.na(group)
counts <- counts[, keep_cells, drop = FALSE]
group <- group[keep_cells]

gene_order <- read.table(
  gene_order_file,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("gene", "chr", "start", "end")
)
gene_order <- gene_order[!duplicated(gene_order$gene), ]
common_genes <- intersect(rownames(counts), gene_order$gene)
if (length(common_genes) < 5000) {
  stop("Too few genes overlap count matrix and gene-order file: ", length(common_genes))
}
counts <- counts[common_genes, , drop = FALSE]

annotations_file <- file.path(input_dir, "GTE009_infercnv_annotations.tsv")
write.table(
  data.frame(cell = colnames(counts), group = group),
  file = annotations_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

summary_file <- file.path(input_dir, "GTE009_infercnv_input_summary.tsv")
summary_table <- data.frame(
  metric = c("cells", "genes", "reference_groups", "tumor_groups", "threads"),
  value = c(
    ncol(counts),
    nrow(counts),
    paste(sort(unique(group[grepl("^Ref_", group)])), collapse = ","),
    paste(sort(unique(group[!grepl("^Ref_", group)])), collapse = ","),
    threads
  )
)
write.table(summary_table, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

message("inferCNV input summary:")
print(table(group))
message("Genes retained: ", nrow(counts))
message("Cells retained: ", ncol(counts))

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = annotations_file,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = sort(unique(group[grepl("^Ref_", group)]))
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = out_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE,
  num_threads = threads,
  plot_steps = FALSE,
  no_prelim_plot = TRUE,
  no_plot = FALSE
)

saveRDS(infercnv_obj, file.path(out_dir, "GTE009_infercnv_object.rds"))
message("inferCNV completed: ", out_dir)