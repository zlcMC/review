#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source("workspace_paths.R")

args <- commandArgs(trailingOnly = TRUE)
infercnv_dir <- if (length(args) >= 1) args[[1]] else output_path("fig2_infercnv_strict", "GTE009")
gene_order_file <- if (length(args) >= 2) args[[2]] else output_path("reference", "gencode_hg38_gene_order.tsv")
out_dir <- if (length(args) >= 3) args[[3]] else output_path("fig2_infercnv_strict", "GTE009_downstream")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

find_first_existing <- function(paths, label) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    stop("Missing ", label, ". Tried:\n  ", paste(paths, collapse = "\n  "))
  }
  existing[[1]]
}

find_optional <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  existing[[1]]
}

read_observations_text <- function(obs_file) {
  message("Reading observations matrix: ", obs_file)
  obs_dt <- fread(obs_file, data.table = FALSE, check.names = FALSE)
  first_col <- names(obs_dt)[[1]]
  if (!first_col %in% colnames(obs_dt)[-1]) {
    gene_ids <- obs_dt[[1]]
    obs_dt[[1]] <- NULL
  } else {
    stop("Could not detect gene-id column in observations matrix: ", obs_file)
  }
  obs_dt[] <- lapply(obs_dt, as.numeric)
  obs_mat <- as.matrix(obs_dt)
  rownames(obs_mat) <- make.unique(as.character(gene_ids))
  obs_mat
}

read_observations_object <- function(object_file, annotations) {
  message("Reading inferCNV object: ", object_file)
  infercnv_obj <- readRDS(object_file)
  if (!isS4(infercnv_obj) || !"expr.data" %in% slotNames(infercnv_obj)) {
    stop("Object does not look like an inferCNV object with an expr.data slot: ", object_file)
  }
  expr_mat <- infercnv_obj@expr.data
  obs_cells <- annotations$cell[annotations$group %in% c("Subclone_1", "Subclone_2")]
  obs_cells <- intersect(obs_cells, colnames(expr_mat))
  if (length(obs_cells) < 20) {
    stop("Too few annotated observation cells found in inferCNV object: ", length(obs_cells))
  }
  expr_mat[, obs_cells, drop = FALSE]
}

annotation_file <- find_first_existing(
  c(
    file.path(infercnv_dir, "input", "GTE009_infercnv_annotations.tsv"),
    file.path(infercnv_dir, "GTE009_infercnv_annotations.tsv")
  ),
  "inferCNV annotation file"
)
if (!file.exists(gene_order_file)) {
  stop("Missing gene-order file: ", gene_order_file)
}

annotations <- fread(annotation_file, header = FALSE, data.table = FALSE)
colnames(annotations) <- c("cell", "group")
annotations$group <- as.character(annotations$group)

obs_file <- find_optional(file.path(infercnv_dir, c(
  "infercnv.observations.txt",
  "infercnv.observations.txt.gz",
  "infercnv.observations.tsv",
  "infercnv.observations.tsv.gz"
)))
if (!is.na(obs_file)) {
  obs_mat <- read_observations_text(obs_file)
} else {
  object_file <- find_first_existing(
    file.path(infercnv_dir, c(
      "GTE009_infercnv_object.rds",
      "run.final.infercnv_obj",
      "22_denoise.leiden.NF_NA.SD_1.5.NL_FALSE.infercnv_obj"
    )),
    "inferCNV final object"
  )
  obs_mat <- read_observations_object(object_file, annotations)
}

annotations <- annotations[annotations$cell %in% colnames(obs_mat), , drop = FALSE]

group1 <- "Subclone_1"
group2 <- "Subclone_2"
cells1 <- annotations$cell[annotations$group == group1]
cells2 <- annotations$cell[annotations$group == group2]
if (length(cells1) < 10 || length(cells2) < 10) {
  stop("Too few cells for downstream comparison: ", group1, "=", length(cells1), ", ", group2, "=", length(cells2))
}

message("Cells: ", group1, "=", length(cells1), "; ", group2, "=", length(cells2))
mat1 <- obs_mat[, cells1, drop = FALSE]
mat2 <- obs_mat[, cells2, drop = FALSE]

row_var <- function(mat) {
  n <- ncol(mat)
  if (n <= 1) {
    return(rep(NA_real_, nrow(mat)))
  }
  means <- rowMeans(mat)
  (rowSums(mat * mat) - n * means * means) / (n - 1)
}

mean1 <- rowMeans(mat1)
mean2 <- rowMeans(mat2)
var1 <- row_var(mat1)
var2 <- row_var(mat2)
n1 <- ncol(mat1)
n2 <- ncol(mat2)
delta <- mean1 - mean2
se <- sqrt(var1 / n1 + var2 / n2)
welch_df <- (var1 / n1 + var2 / n2)^2 / ((var1 / n1)^2 / (n1 - 1) + (var2 / n2)^2 / (n2 - 1))
t_stat <- delta / se
p_value <- 2 * pt(abs(t_stat), df = welch_df, lower.tail = FALSE)
p_value[!is.finite(p_value)] <- NA_real_

stats <- data.frame(
  gene = rownames(obs_mat),
  mean_subclone1 = mean1,
  mean_subclone2 = mean2,
  delta_subclone1_minus_subclone2 = delta,
  t_stat = t_stat,
  p_value = p_value,
  fdr = p.adjust(p_value, method = "BH"),
  stringsAsFactors = FALSE
)

gene_order <- fread(gene_order_file, header = FALSE, data.table = FALSE)
colnames(gene_order) <- c("gene", "chr", "start", "end")
gene_order <- gene_order[!duplicated(gene_order$gene), ]
stats <- merge(stats, gene_order, by = "gene", all.x = TRUE, sort = FALSE)
stats <- stats[order(stats$p_value), ]
write.csv(stats, file.path(out_dir, "GTE009_inferCNV_genelevel_subclone_stats.csv"), row.names = FALSE)

dyn <- stats[stats$gene == "DYNC2H1", , drop = FALSE]
write.csv(dyn, file.path(out_dir, "GTE009_inferCNV_DYNC2H1_stats.csv"), row.names = FALSE)
if (nrow(dyn) == 0) {
  warning("DYNC2H1 was not found in inferCNV observations. Volcano will still be generated.")
} else {
  message("DYNC2H1 delta: ", signif(dyn$delta_subclone1_minus_subclone2, 4), "; p=", signif(dyn$p_value, 4))
}

plot_df <- stats
plot_df$neglog10_fdr <- -log10(pmax(plot_df$fdr, .Machine$double.xmin))
plot_df$is_dync2h1 <- plot_df$gene == "DYNC2H1"
plot_df$is_chr11 <- plot_df$chr == "chr11"
plot_df$chromosome_set <- ifelse(plot_df$is_chr11, "chr11", "Other chromosomes")
plot_df$label <- ifelse(plot_df$is_dync2h1, "DYNC2H1", NA_character_)

p_volcano <- ggplot(plot_df, aes(delta_subclone1_minus_subclone2, neglog10_fdr)) +
  geom_point(aes(color = chromosome_set), size = 0.8, alpha = 0.45) +
  geom_point(data = subset(plot_df, is_dync2h1), color = "#D7191C", size = 2.8) +
  geom_text(data = subset(plot_df, is_dync2h1), aes(label = label), vjust = -0.8, color = "#D7191C", size = 3.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey45") +
  scale_color_manual(values = c(`chr11` = "#2C7BB6", `Other chromosomes` = "grey70"), name = NULL) +
  labs(
    title = "GTE009 inferCNV gene-level difference",
    subtitle = "Positive delta indicates higher inferred CNV signal in Subclone 1",
    x = "Mean inferCNV signal: Subclone 1 - Subclone 2",
    y = "-log10(FDR)"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave(file.path(out_dir, "Fig2D_inferCNV_DYNC2H1_genelevel_volcano.pdf"), p_volcano, width = 7.2, height = 5.2)
ggsave(file.path(out_dir, "Fig2D_inferCNV_DYNC2H1_genelevel_volcano.png"), p_volcano, width = 7.2, height = 5.2, dpi = 260)

chr11_genes <- gene_order$gene[gene_order$chr == "chr11"]
chr11_genes <- intersect(chr11_genes, rownames(obs_mat))
if (length(chr11_genes) < 20) {
  stop("Too few chr11 genes found in inferCNV observations: ", length(chr11_genes))
}

chr11_stats <- stats[match(chr11_genes, stats$gene), ]
chr11_stats <- chr11_stats[order(chr11_stats$start), ]
chr11_stats$position_mb <- (chr11_stats$start + chr11_stats$end) / 2 / 1e6
chr11_stats$gene_index <- seq_len(nrow(chr11_stats))
chr11_long <- rbind(
  data.frame(gene = chr11_stats$gene, gene_index = chr11_stats$gene_index, position_mb = chr11_stats$position_mb, group = group1, mean_cnv = chr11_stats$mean_subclone1),
  data.frame(gene = chr11_stats$gene, gene_index = chr11_stats$gene_index, position_mb = chr11_stats$position_mb, group = group2, mean_cnv = chr11_stats$mean_subclone2)
)
chr11_long$group <- factor(chr11_long$group, levels = c(group1, group2))

axis_mbs <- seq(0, floor(max(chr11_stats$position_mb, na.rm = TRUE) / 25) * 25, by = 25)
axis_idx <- vapply(axis_mbs, function(mb) {
  chr11_stats$gene_index[which.min(abs(chr11_stats$position_mb - mb))]
}, numeric(1))
dync_chr11 <- subset(chr11_stats, gene == "DYNC2H1")

midpoint <- 1
p_chr11_mean <- ggplot(chr11_long, aes(gene_index, group, fill = mean_cnv)) +
  geom_tile(width = 1, height = 0.86) +
  geom_vline(data = dync_chr11, aes(xintercept = gene_index), color = "#D7191C", linewidth = 0.55) +
  scale_x_continuous(breaks = axis_idx, labels = axis_mbs, expand = c(0.005, 0.005)) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = midpoint, name = "inferCNV") +
  labs(
    title = "GTE009 chr11 inferCNV signal by subclone",
    subtitle = "Vertical red line marks DYNC2H1",
    x = "chr11 position (Mb; genes ordered by genomic coordinate)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))

ggsave(file.path(out_dir, "Fig2E_chr11_inferCNV_group_mean_heatmap.pdf"), p_chr11_mean, width = 8.2, height = 2.8)
ggsave(file.path(out_dir, "Fig2E_chr11_inferCNV_group_mean_heatmap.png"), p_chr11_mean, width = 8.2, height = 2.8, dpi = 260)
write.csv(chr11_stats, file.path(out_dir, "GTE009_chr11_inferCNV_genelevel_stats.csv"), row.names = FALSE)

p_chr11_delta <- ggplot(chr11_stats, aes(position_mb, delta_subclone1_minus_subclone2)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey45") +
  geom_line(color = "grey40", linewidth = 0.55) +
  geom_point(data = dync_chr11, color = "#D7191C", size = 2.6) +
  geom_text(data = dync_chr11, aes(label = gene), vjust = -0.45, color = "#D7191C", size = 3.4) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  labs(
    title = "GTE009 chr11 inferCNV subclone delta profile",
    subtitle = "Positive values indicate higher inferred CNV signal in Subclone 1",
    x = "chr11 position (Mb)",
    y = "Subclone 1 - Subclone 2"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "Fig2E_chr11_inferCNV_delta_profile.pdf"), p_chr11_delta, width = 8.2, height = 3.2)
ggsave(file.path(out_dir, "Fig2E_chr11_inferCNV_delta_profile.png"), p_chr11_delta, width = 8.2, height = 3.2, dpi = 260)

heatmap_grob <- NULL
if (requireNamespace("pheatmap", quietly = TRUE)) {
  set.seed(42)
  max_cells_per_group <- 500
  selected_cells <- c(
    sample(cells1, min(max_cells_per_group, length(cells1))),
    sample(cells2, min(max_cells_per_group, length(cells2)))
  )
  max_chr11_genes <- 700
  selected_chr11 <- chr11_stats$gene
  if (length(selected_chr11) > max_chr11_genes) {
    sampled_chr11 <- selected_chr11[round(seq(1, length(selected_chr11), length.out = max_chr11_genes))]
    selected_set <- unique(c(sampled_chr11, "DYNC2H1"))
    selected_chr11 <- chr11_stats$gene[chr11_stats$gene %in% selected_set]
  }

  heat_mat <- obs_mat[selected_chr11, selected_cells, drop = FALSE]
  heat_mat <- heat_mat[match(selected_chr11, rownames(heat_mat)), , drop = FALSE]
  col_annot <- data.frame(group = annotations$group[match(colnames(heat_mat), annotations$cell)])
  rownames(col_annot) <- colnames(heat_mat)
  heat_mat <- heat_mat[, order(col_annot$group), drop = FALSE]
  col_annot <- col_annot[colnames(heat_mat), , drop = FALSE]
  breaks <- seq(quantile(heat_mat, 0.02, na.rm = TRUE), quantile(heat_mat, 0.98, na.rm = TRUE), length.out = 101)
  colors <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100)
  heatmap_grob <- pheatmap::pheatmap(
    heat_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_col = col_annot,
    color = colors,
    breaks = breaks,
    main = "GTE009 chr11 inferCNV sampled cells",
    silent = TRUE
  )$gtable
  pdf(file.path(out_dir, "Fig2F_chr11_inferCNV_sampled_cell_heatmap.pdf"), width = 10, height = 7)
  grid::grid.newpage()
  grid::grid.draw(heatmap_grob)
  dev.off()
  png(file.path(out_dir, "Fig2F_chr11_inferCNV_sampled_cell_heatmap.png"), width = 2600, height = 1800, res = 260)
  grid::grid.newpage()
  grid::grid.draw(heatmap_grob)
  dev.off()
}

if (!is.null(heatmap_grob)) {
  panel_tag_theme <- theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.98)
  )
  p_volcano_panel <- p_volcano + labs(tag = "D") + panel_tag_theme
  p_chr11_mean_panel <- p_chr11_mean + labs(tag = "E") + panel_tag_theme
  p_chr11_delta_panel <- p_chr11_delta + labs(tag = "E'") + panel_tag_theme

  draw_combined_panel <- function() {
    grid::grid.newpage()
    layout <- grid::grid.layout(
      nrow = 3,
      ncol = 2,
      heights = grid::unit(c(0.95, 0.68, 1.35), "null"),
      widths = grid::unit(c(1, 1), "null")
    )
    grid::pushViewport(grid::viewport(layout = layout))
    print(p_volcano_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_chr11_mean_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(p_chr11_delta_panel, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 1:2))
    grid::grid.draw(heatmap_grob)
    grid::grid.text("F", x = grid::unit(0.01, "npc"), y = grid::unit(0.98, "npc"), just = c("left", "top"), gp = grid::gpar(fontface = "bold", fontsize = 16))
    grid::popViewport()
    grid::popViewport()
  }

  pdf(file.path(out_dir, "Fig2D_F_inferCNV_combined_panel.pdf"), width = 13.5, height = 12.2)
  draw_combined_panel()
  dev.off()
  png(file.path(out_dir, "Fig2D_F_inferCNV_combined_panel.png"), width = 3510, height = 3172, res = 260)
  draw_combined_panel()
  dev.off()
}

message("Downstream inferCNV figures written to: ", out_dir)