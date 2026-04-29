#!/usr/bin/env Rscript
# v2 of projectmd/fig02_infercnv_genelevel_panels.R
# 与原脚本逻辑/默认参数完全一致；切换到 v2 路径 helper。
# 默认从 output/v2 读 inferCNV 主对象，缺失时回退到 v1 (output/fig2_infercnv_strict/GTE009/)。
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
this_script <- sub('^--file=', '',
                   grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(this_script), "helpers", "helpers.R"))
init_workspace()

args <- commandArgs(trailingOnly = TRUE)
infercnv_dir    <- if (length(args) >= 1) args[[1]] else dep_file("fig2_infercnv_strict", "GTE009")
gene_order_file <- if (length(args) >= 2) args[[2]] else dep_file("reference", "gencode_hg38_gene_order.tsv")
out_dir         <- if (length(args) >= 3) args[[3]] else v2_dir("fig2_infercnv_strict", "GTE009_downstream")

find_first <- function(paths, label) {
  ok <- paths[file.exists(paths)]
  if (!length(ok)) stop("Missing ", label, ":\n  ", paste(paths, collapse = "\n  "))
  ok[[1]]
}
find_optional <- function(paths) { ok <- paths[file.exists(paths)]; if (length(ok)) ok[[1]] else NA_character_ }

read_obs_text <- function(f) {
  message("Reading observations: ", f)
  d <- fread(f, data.table = FALSE, check.names = FALSE)
  g <- d[[1]]; d[[1]] <- NULL
  d[] <- lapply(d, as.numeric)
  m <- as.matrix(d); rownames(m) <- make.unique(as.character(g)); m
}
read_obs_obj <- function(f, ann) {
  message("Reading inferCNV object: ", f)
  o <- readRDS(f)
  stopifnot(isS4(o), "expr.data" %in% slotNames(o))
  e <- o@expr.data
  cells <- intersect(ann$cell[ann$group %in% c("Subclone_1", "Subclone_2")], colnames(e))
  stopifnot(length(cells) >= 20)
  e[, cells, drop = FALSE]
}

ann_file <- find_first(c(file.path(infercnv_dir, "input", "GTE009_infercnv_annotations.tsv"),
                          file.path(infercnv_dir, "GTE009_infercnv_annotations.tsv")),
                        "annotation file")
stopifnot(file.exists(gene_order_file))

ann <- fread(ann_file, header = FALSE, data.table = FALSE)
colnames(ann) <- c("cell", "group"); ann$group <- as.character(ann$group)

obs_file <- find_optional(file.path(infercnv_dir,
              c("infercnv.observations.txt", "infercnv.observations.txt.gz",
                "infercnv.observations.tsv", "infercnv.observations.tsv.gz")))
obs_mat <- if (!is.na(obs_file)) read_obs_text(obs_file) else
  read_obs_obj(find_first(file.path(infercnv_dir,
              c("GTE009_infercnv_object.rds", "run.final.infercnv_obj",
                "22_denoise.leiden.NF_NA.SD_1.5.NL_FALSE.infercnv_obj")),
              "inferCNV final object"), ann)

ann <- ann[ann$cell %in% colnames(obs_mat), , drop = FALSE]
g1 <- "Subclone_1"; g2 <- "Subclone_2"
c1 <- ann$cell[ann$group == g1]; c2 <- ann$cell[ann$group == g2]
stopifnot(length(c1) >= 10, length(c2) >= 10)
message("Cells: ", g1, "=", length(c1), "; ", g2, "=", length(c2))

m1 <- obs_mat[, c1, drop = FALSE]; m2 <- obs_mat[, c2, drop = FALSE]
row_var <- function(M) {
  n <- ncol(M); if (n <= 1) return(rep(NA_real_, nrow(M)))
  mu <- rowMeans(M); (rowSums(M * M) - n * mu * mu) / (n - 1)
}
mu1 <- rowMeans(m1); mu2 <- rowMeans(m2); v1 <- row_var(m1); v2 <- row_var(m2)
n1 <- ncol(m1); n2 <- ncol(m2)
delta <- mu1 - mu2; se <- sqrt(v1 / n1 + v2 / n2)
df    <- (v1 / n1 + v2 / n2)^2 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
tstat <- delta / se
pval  <- 2 * pt(abs(tstat), df = df, lower.tail = FALSE); pval[!is.finite(pval)] <- NA_real_

stats <- data.frame(gene = rownames(obs_mat),
                    mean_subclone1 = mu1, mean_subclone2 = mu2,
                    delta_subclone1_minus_subclone2 = delta,
                    t_stat = tstat, p_value = pval,
                    fdr = p.adjust(pval, method = "BH"),
                    stringsAsFactors = FALSE)
go <- fread(gene_order_file, header = FALSE, data.table = FALSE)
colnames(go) <- c("gene", "chr", "start", "end")
go <- go[!duplicated(go$gene), ]
stats <- merge(stats, go, by = "gene", all.x = TRUE, sort = FALSE)
stats <- stats[order(stats$p_value), ]
write.csv(stats, file.path(out_dir, "GTE009_inferCNV_genelevel_subclone_stats.csv"), row.names = FALSE)

dyn <- stats[stats$gene == "DYNC2H1", , drop = FALSE]
write.csv(dyn, file.path(out_dir, "GTE009_inferCNV_DYNC2H1_stats.csv"), row.names = FALSE)
if (nrow(dyn)) message("DYNC2H1 delta=", signif(dyn$delta_subclone1_minus_subclone2, 4),
                       " p=", signif(dyn$p_value, 4))

pdf_png <- function(p, base, w, h, dpi = 260) {
  ggsave(paste0(base, ".pdf"), p, width = w, height = h)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = dpi)
}

plot_df <- transform(stats,
  neglog10_fdr   = -log10(pmax(fdr, .Machine$double.xmin)),
  is_dync2h1     = gene == "DYNC2H1",
  chromosome_set = ifelse(chr == "chr11", "chr11", "Other chromosomes"),
  label          = ifelse(gene == "DYNC2H1", "DYNC2H1", NA_character_))

p_volcano <- ggplot(plot_df, aes(delta_subclone1_minus_subclone2, neglog10_fdr)) +
  geom_point(aes(color = chromosome_set), size = 0.8, alpha = 0.45) +
  geom_point(data = subset(plot_df, is_dync2h1), color = "#D7191C", size = 2.8) +
  geom_text(data = subset(plot_df, is_dync2h1), aes(label = label),
            vjust = -0.8, color = "#D7191C", size = 3.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey45") +
  scale_color_manual(values = c(`chr11` = "#2C7BB6", `Other chromosomes` = "grey70"), name = NULL) +
  labs(title = "GTE009 inferCNV gene-level difference",
       subtitle = "Positive delta indicates higher inferred CNV signal in Subclone 1",
       x = "Mean inferCNV signal: Subclone 1 - Subclone 2", y = "-log10(FDR)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")
pdf_png(p_volcano, file.path(out_dir, "Fig2D_inferCNV_DYNC2H1_genelevel_volcano"), 7.2, 5.2)

chr11_genes <- intersect(go$gene[go$chr == "chr11"], rownames(obs_mat))
stopifnot(length(chr11_genes) >= 20)
chr11_stats <- stats[match(chr11_genes, stats$gene), ]
chr11_stats <- chr11_stats[order(chr11_stats$start), ]
chr11_stats$position_mb <- (chr11_stats$start + chr11_stats$end) / 2 / 1e6
chr11_stats$gene_index  <- seq_len(nrow(chr11_stats))
write.csv(chr11_stats, file.path(out_dir, "GTE009_chr11_inferCNV_genelevel_stats.csv"), row.names = FALSE)

chr11_long <- rbind(
  data.frame(gene = chr11_stats$gene, gene_index = chr11_stats$gene_index,
             position_mb = chr11_stats$position_mb, group = g1, mean_cnv = chr11_stats$mean_subclone1),
  data.frame(gene = chr11_stats$gene, gene_index = chr11_stats$gene_index,
             position_mb = chr11_stats$position_mb, group = g2, mean_cnv = chr11_stats$mean_subclone2))
chr11_long$group <- factor(chr11_long$group, levels = c(g1, g2))

axis_mbs <- seq(0, floor(max(chr11_stats$position_mb, na.rm = TRUE) / 25) * 25, by = 25)
axis_idx <- vapply(axis_mbs,
                   function(mb) chr11_stats$gene_index[which.min(abs(chr11_stats$position_mb - mb))],
                   numeric(1))
dync <- subset(chr11_stats, gene == "DYNC2H1")

p_chr11_mean <- ggplot(chr11_long, aes(gene_index, group, fill = mean_cnv)) +
  geom_tile(width = 1, height = 0.86) +
  geom_vline(data = dync, aes(xintercept = gene_index), color = "#D7191C", linewidth = 0.55) +
  scale_x_continuous(breaks = axis_idx, labels = axis_mbs, expand = c(0.005, 0.005)) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 1, name = "inferCNV") +
  labs(title = "GTE009 chr11 inferCNV signal by subclone",
       subtitle = "Vertical red line marks DYNC2H1",
       x = "chr11 position (Mb; genes ordered by genomic coordinate)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))
pdf_png(p_chr11_mean, file.path(out_dir, "Fig2E_chr11_inferCNV_group_mean_heatmap"), 8.2, 2.8)

p_chr11_delta <- ggplot(chr11_stats, aes(position_mb, delta_subclone1_minus_subclone2)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey45") +
  geom_line(color = "grey40", linewidth = 0.55) +
  geom_point(data = dync, color = "#D7191C", size = 2.6) +
  geom_text(data = dync, aes(label = gene), vjust = -0.45, color = "#D7191C", size = 3.4) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  labs(title = "GTE009 chr11 inferCNV subclone delta profile",
       subtitle = "Positive values indicate higher inferred CNV signal in Subclone 1",
       x = "chr11 position (Mb)", y = "Subclone 1 - Subclone 2") +
  theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold"))
pdf_png(p_chr11_delta, file.path(out_dir, "Fig2E_chr11_inferCNV_delta_profile"), 8.2, 3.2)

heatmap_grob <- NULL
if (requireNamespace("pheatmap", quietly = TRUE)) {
  set.seed(42)
  sel_cells <- c(sample(c1, min(500, length(c1))), sample(c2, min(500, length(c2))))
  sel_chr11 <- chr11_stats$gene
  if (length(sel_chr11) > 700) {
    samp <- sel_chr11[round(seq(1, length(sel_chr11), length.out = 700))]
    sel_chr11 <- chr11_stats$gene[chr11_stats$gene %in% unique(c(samp, "DYNC2H1"))]
  }
  hm <- obs_mat[sel_chr11, sel_cells, drop = FALSE]
  hm <- hm[match(sel_chr11, rownames(hm)), , drop = FALSE]
  ca <- data.frame(group = ann$group[match(colnames(hm), ann$cell)])
  rownames(ca) <- colnames(hm); hm <- hm[, order(ca$group), drop = FALSE]
  ca <- ca[colnames(hm), , drop = FALSE]
  brks <- seq(quantile(hm, 0.02, na.rm = TRUE), quantile(hm, 0.98, na.rm = TRUE), length.out = 101)
  cols <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100)
  heatmap_grob <- pheatmap::pheatmap(hm, cluster_rows = FALSE, cluster_cols = FALSE,
    show_rownames = FALSE, show_colnames = FALSE, annotation_col = ca,
    color = cols, breaks = brks, main = "GTE009 chr11 inferCNV sampled cells",
    silent = TRUE)$gtable
  pdf(file.path(out_dir, "Fig2F_chr11_inferCNV_sampled_cell_heatmap.pdf"), width = 10, height = 7)
  grid::grid.newpage(); grid::grid.draw(heatmap_grob); dev.off()
  png(file.path(out_dir, "Fig2F_chr11_inferCNV_sampled_cell_heatmap.png"),
      width = 2600, height = 1800, res = 260)
  grid::grid.newpage(); grid::grid.draw(heatmap_grob); dev.off()
}

if (!is.null(heatmap_grob)) {
  tag_theme <- theme(plot.tag = element_text(face = "bold", size = 16),
                     plot.tag.position = c(0.01, 0.98))
  pV <- p_volcano     + labs(tag = "D")  + tag_theme
  pM <- p_chr11_mean  + labs(tag = "E")  + tag_theme
  pD <- p_chr11_delta + labs(tag = "E'") + tag_theme

  draw_combined <- function() {
    grid::grid.newpage()
    lay <- grid::grid.layout(nrow = 3, ncol = 2,
      heights = grid::unit(c(0.95, 0.68, 1.35), "null"),
      widths  = grid::unit(c(1, 1), "null"))
    grid::pushViewport(grid::viewport(layout = lay))
    print(pV, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(pM, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(pD, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 1:2))
    grid::grid.draw(heatmap_grob)
    grid::grid.text("F", x = grid::unit(0.01, "npc"), y = grid::unit(0.98, "npc"),
                    just = c("left", "top"), gp = grid::gpar(fontface = "bold", fontsize = 16))
    grid::popViewport(); grid::popViewport()
  }
  pdf(file.path(out_dir, "Fig2D_F_inferCNV_combined_panel.pdf"), width = 13.5, height = 12.2)
  draw_combined(); dev.off()
  png(file.path(out_dir, "Fig2D_F_inferCNV_combined_panel.png"),
      width = 3510, height = 3172, res = 260)
  draw_combined(); dev.off()
}
