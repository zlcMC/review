#!/usr/bin/env Rscript
# v2 of fig05_cellchat_readable_figures.R — reader-friendly heatmaps/dotplots.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(CellChat); library(ggplot2); library(patchwork)
})

out    <- v2_dir('fig5_readable')
cc     <- cellchat_load()
g      <- cellchat_groups(cc)
groups <- g$all
ct_pal <- cellchat_palette(groups)

message('Groups: ', paste(groups, collapse = ', '))
message('Malignant-like: ', paste(g$mal, collapse = ', '))
message('Microenvironment: ', paste(g$micro, collapse = ', '))

comm_all <- subsetCommunication(cc)
safe_prob_limits <- range(comm_all$prob, na.rm = TRUE)
prob_cols <- c('#315CA8','#37A7B8','#86D19E','#F6D96B','#E86A45','#B3133B')

theme_readable <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y  = element_text(size = base_size),
          axis.title   = element_text(size = base_size + 1),
          plot.title   = element_text(size = base_size + 3, face = 'bold', hjust = 0),
          legend.title = element_text(size = base_size),
          legend.text  = element_text(size = base_size - 1),
          plot.margin  = margin(8, 12, 8, 12))
}

# Helper: return top-N (pathway × interaction) by max prob.
top_interactions <- function(comm, top_n) {
  rk <- aggregate(prob ~ pathway_name + interaction_name, data = comm, FUN = max)
  rk <- rk[order(-rk$prob, rk$pathway_name, rk$interaction_name), , drop = FALSE]
  rk <- head(rk, top_n)
  rk$interaction_label <- paste(rk$pathway_name, rk$interaction_name, sep = ': ')
  rk
}

plot_direction_bubble <- function(comm, sources, targets, direction, title, top_n = 45) {
  st_levels <- as.vector(outer(sources, targets, paste, sep = ' -> '))
  comm$source_target <- paste(comm$source, comm$target, sep = ' -> ')
  comm <- comm[comm$source_target %in% st_levels, , drop = FALSE]
  if (!nrow(comm)) { warning('No rows for ', direction); return(invisible()) }

  rk <- top_interactions(comm, top_n)
  comm <- merge(comm, rk[, c('pathway_name','interaction_name','interaction_label')],
                by = c('pathway_name','interaction_name'))
  comm$source_target     <- factor(comm$source_target,     levels = st_levels)
  comm$interaction_label <- factor(comm$interaction_label, levels = rev(rk$interaction_label))
  comm$p_group <- ifelse(comm$pval < 0.01, 'p < 0.01', '0.01 <= p < 0.05')

  write.csv(rk, file.path(out, paste0('Fig5B_', direction, '_top_interactions.csv')),
            row.names = FALSE)

  w <- max(12, length(st_levels) * 0.22 + 4)
  h <- max(8,  length(unique(comm$interaction_label)) * 0.18 + 2)
  ggsave(file.path(out, paste0('Fig5B_bubble_', direction, '_top45.pdf')),
         ggplot(comm, aes(x = source_target, y = interaction_label)) +
           geom_point(aes(color = prob, size = p_group), alpha = 0.9) +
           scale_color_gradientn(colors = prob_cols, limits = safe_prob_limits,
                                 name = 'Communication\nprobability') +
           scale_size_manual(values = c('0.01 <= p < 0.05' = 1.2, 'p < 0.01' = 2.4),
                             name = 'Significance') +
           labs(title = title, x = 'Source -> target cell group',
                y = 'Pathway: ligand-receptor pair') +
           theme_readable(7) +
           theme(panel.grid.major = element_line(color = 'grey88', linewidth = 0.25),
                 axis.text.y      = element_text(size = 5.7)),
         width = w, height = h, limitsize = FALSE)
}

comm_mal_to_micro <- subsetCommunication(cc, sources.use = g$mal,   targets.use = g$micro)
comm_micro_to_mal <- subsetCommunication(cc, sources.use = g$micro, targets.use = g$mal)
plot_direction_bubble(comm_mal_to_micro, g$mal,   g$micro, 'mal_to_micro',
                      'Top ligand-receptor signals: malignant-like -> microenvironment')
plot_direction_bubble(comm_micro_to_mal, g$micro, g$mal,   'micro_to_mal',
                      'Top ligand-receptor signals: microenvironment -> malignant-like')

plot_split_bubble <- function(comm, fixed_col, fixed_vals, vary_col,
                              direction, title_prefix, x_label, top_n = 30) {
  for (v in fixed_vals) {
    fc <- comm[comm[[fixed_col]] == v, , drop = FALSE]
    if (!nrow(fc)) next
    rk <- top_interactions(fc, top_n)
    fc <- merge(fc, rk[, c('pathway_name','interaction_name','interaction_label')],
                by = c('pathway_name','interaction_name'))
    fc$cell_group        <- factor(fc[[vary_col]], levels = groups)
    fc$interaction_label <- factor(fc$interaction_label, levels = rev(rk$interaction_label))
    fc$p_group           <- ifelse(fc$pval < 0.01, 'p < 0.01', '0.01 <= p < 0.05')

    stem <- paste0('Fig5B_bubble_', direction, '_', v, '_top30')
    write.csv(rk, file.path(out, paste0(stem, '.csv')), row.names = FALSE)
    h <- max(7, length(unique(fc$interaction_label)) * 0.19 + 2)
    ggsave(file.path(out, paste0(stem, '.pdf')),
           ggplot(fc, aes(x = cell_group, y = interaction_label)) +
             geom_point(aes(color = prob, size = p_group), alpha = 0.9) +
             scale_color_gradientn(colors = prob_cols, limits = safe_prob_limits,
                                   name = 'Communication\nprobability') +
             scale_size_manual(values = c('0.01 <= p < 0.05' = 1.3, 'p < 0.01' = 2.8),
                               name = 'Significance') +
             labs(title = paste0(title_prefix, v), x = x_label,
                  y = 'Pathway: ligand-receptor pair') +
             theme_readable(8) +
             theme(panel.grid.major = element_line(color = 'grey88', linewidth = 0.25),
                   axis.text.y      = element_text(size = 6.3),
                   axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1)),
           width = 8.2, height = h, limitsize = FALSE)
  }
}
plot_split_bubble(comm_mal_to_micro, 'target', g$micro, 'source',
                  'mal_to_micro_target', 'Malignant-like -> ',
                  'Source malignant-like cell group')
plot_split_bubble(comm_micro_to_mal, 'source', g$micro, 'target',
                  'micro_to_mal_source', 'Microenvironment source: ',
                  'Target malignant-like cell group')

# Overall crosstalk: heatmaps + top-edge bars.
plot_matrix_heatmap <- function(mat, title, file, w = 8.5, h = 7.5) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c('source','target','strength')
  df$source <- factor(df$source, levels = rev(rownames(mat)))
  df$target <- factor(df$target, levels = colnames(mat))
  ggsave(file.path(out, file),
         ggplot(df, aes(x = target, y = source, fill = strength)) +
           geom_tile(color = 'white', linewidth = 0.25) +
           scale_fill_gradientn(colors = c('#FFF7EC','#FDD49E','#FC8D59','#D7301F','#7F0000'),
                                name = 'Strength') +
           labs(title = title, x = 'Receiver / target', y = 'Sender / source') +
           coord_fixed() + theme_readable(7) +
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.border = element_rect(color = 'grey50', fill = NA, linewidth = 0.3)),
         width = w, height = h, limitsize = FALSE)
}

plot_top_edges <- function(mat, title, file, top_n = 40) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c('source','target','strength')
  df <- df[df$strength > 0, , drop = FALSE]
  df <- head(df[order(-df$strength), , drop = FALSE], top_n)
  df$edge <- factor(paste(df$source, df$target, sep = ' -> '),
                    levels = rev(paste(df$source, df$target, sep = ' -> ')))
  write.csv(df, file.path(out, sub('pdf$', 'csv', file)), row.names = FALSE)
  ggsave(file.path(out, file),
         ggplot(df, aes(x = strength, y = edge, fill = source)) +
           geom_col(width = 0.78, show.legend = FALSE) +
           scale_fill_manual(values = ct_pal, na.value = '#999999') +
           labs(title = title, x = 'Aggregated interaction strength',
                y = 'Source -> target') +
           theme_classic(8) +
           theme(plot.title  = element_text(size = 11, face = 'bold'),
                 axis.text.y = element_text(size = 6.5),
                 plot.margin = margin(8, 12, 8, 12)),
         width = 8.3, height = 9.2, limitsize = FALSE)
}

plot_matrix_heatmap(cc@net$weight, 'Overall CellChat interaction strength',
                    'Fig5B_crosstalk_strength_heatmap_readable.pdf')
plot_matrix_heatmap(cc@net$count, 'Overall CellChat interaction count',
                    'Fig5B_crosstalk_count_heatmap_readable.pdf')
plot_top_edges(cc@net$weight, 'Top source-target pairs by interaction strength',
               'Fig5B_crosstalk_top40_edges.pdf')

# Pathway-level heatmaps.
plot_pathway_heatmap <- function(pathway, individual = TRUE) {
  mat <- cc@netP$prob[, , pathway]
  df  <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c('source','target','prob')
  df$source <- factor(df$source, levels = rev(rownames(mat)))
  df$target <- factor(df$target, levels = colnames(mat))
  p <- ggplot(df, aes(x = target, y = source, fill = prob)) +
    geom_tile(color = 'white', linewidth = 0.25) +
    scale_fill_gradientn(colors = prob_cols, name = 'Prob.') +
    labs(title = paste0(pathway, ' pathway'),
         x = 'Receiver / target', y = 'Sender / source') +
    coord_fixed() +
    theme_readable(if (individual) 7 else 5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(size = if (individual) 11 else 7, face = 'bold'),
          legend.position = if (individual) 'right' else 'none')
  if (individual) {
    ggsave(file.path(out, paste0('Fig5C_pathway_heatmap_', pathway, '.pdf')),
           p, width = 8.5, height = 7.5, limitsize = FALSE)
  }
  p
}

focus_pw <- intersect(c('MK','PTN','MIF','SPP1','GALECTIN','TGFb','VEGF'), cc@netP$pathways)
invisible(lapply(focus_pw, plot_pathway_heatmap, individual = TRUE))
panel <- wrap_plots(lapply(focus_pw, plot_pathway_heatmap, individual = FALSE), ncol = 4)
ggsave(file.path(out, 'Fig5C_pathway_heatmap_panel_readable.pdf'),
       panel, width = 15, height = 8.5, limitsize = FALSE)

# NMF pattern summaries.
plot_pattern_signaling <- function(direction, file_prefix, top_n = 10) {
  pat <- cc@netP$pattern[[direction]]$pattern
  if (is.null(pat) || is.null(pat$signaling) || is.null(pat$cell)) {
    warning('No pattern data for ', direction); return(invisible())
  }
  sig <- pat$signaling
  sig$Pattern   <- as.character(sig$Pattern)
  sig$Signaling <- as.character(sig$Signaling)
  idx <- unlist(lapply(split(seq_len(nrow(sig)), sig$Pattern), function(i) {
    head(i[order(-sig$Contribution[i])], top_n)
  }), use.names = FALSE)
  sig_top <- sig[idx, , drop = FALSE]
  sig_top <- sig_top[order(sig_top$Pattern, -sig_top$Contribution), , drop = FALSE]
  sig_top$Signaling <- factor(sig_top$Signaling, levels = rev(unique(sig_top$Signaling)))
  write.csv(sig_top, file.path(out, paste0(file_prefix, '_top_signaling.csv')), row.names = FALSE)

  ggsave(file.path(out, paste0(file_prefix, '_top_signaling_dotplot.pdf')),
         ggplot(sig_top, aes(x = Pattern, y = Signaling)) +
           geom_point(aes(size = Contribution, color = Contribution), alpha = 0.92) +
           scale_color_gradientn(colors = prob_cols, name = 'Contribution') +
           scale_size(range = c(1.2, 4.2), name = 'Contribution') +
           labs(title = paste0(direction, ' patterns: top signaling pathways'),
                x = NULL, y = 'Signaling pathway') +
           theme_classic(8) +
           theme(plot.title  = element_text(size = 11, face = 'bold'),
                 axis.text.y = element_text(size = 6.5),
                 panel.grid.major.y = element_line(color = 'grey88', linewidth = 0.25)),
         width = 7.5, height = max(6, length(unique(sig_top$Signaling)) * 0.17 + 2),
         limitsize = FALSE)

  sig_full <- sig
  dom_idx <- unlist(lapply(split(seq_len(nrow(sig_full)), sig_full$Signaling), function(i) {
    i[which.max(sig_full$Contribution[i])]
  }), use.names = FALSE)
  dom <- sig_full[dom_idx, , drop = FALSE]
  dom <- dom[order(dom$Pattern, -dom$Contribution, dom$Signaling), , drop = FALSE]
  sig_full$Pattern   <- factor(sig_full$Pattern, levels = sort(unique(sig_full$Pattern)))
  sig_full$Signaling <- factor(sig_full$Signaling, levels = rev(dom$Signaling))
  write.csv(sig_full, file.path(out, paste0(file_prefix, '_all_signaling.csv')), row.names = FALSE)

  full_h <- max(11, length(unique(sig_full$Signaling)) * 0.13 + 2)
  p_full <- ggplot(sig_full, aes(x = Pattern, y = Signaling, fill = Contribution)) +
    geom_tile(color = 'white', linewidth = 0.18) +
    scale_fill_gradientn(colors = prob_cols, name = 'Contribution') +
    labs(title = paste0(direction, ' patterns: all signaling pathways'),
         x = NULL, y = 'Signaling pathway') +
    theme_classic(8) +
    theme(plot.title  = element_text(size = 11, face = 'bold'),
          axis.text.y = element_text(size = 5.2),
          panel.border = element_rect(color = 'grey50', fill = NA, linewidth = 0.3))
  ggsave(file.path(out, paste0(file_prefix, '_all_signaling_heatmap.pdf')),
         p_full, width = 6.6, height = full_h, limitsize = FALSE)

  cell <- pat$cell
  cell$Pattern   <- as.character(cell$Pattern)
  cell$CellGroup <- factor(as.character(cell$CellGroup), levels = rev(groups))
  p_cell <- ggplot(cell, aes(x = Pattern, y = CellGroup, fill = Contribution)) +
    geom_tile(color = 'white', linewidth = 0.25) +
    scale_fill_gradientn(colors = prob_cols, name = 'Contribution') +
    labs(title = paste0(direction, ' patterns: cell-group contributions'),
         x = NULL, y = NULL) +
    theme_classic(8) +
    theme(plot.title   = element_text(size = 11, face = 'bold'),
          panel.border = element_rect(color = 'grey50', fill = NA, linewidth = 0.3))
  ggsave(file.path(out, paste0(file_prefix, '_cell_contribution_heatmap.pdf')),
         p_cell, width = 5.8, height = 5.8, limitsize = FALSE)
  ggsave(file.path(out, paste0(file_prefix, '_readable_overview.pdf')),
         p_cell + p_full + plot_layout(widths = c(0.95, 1.35)),
         width = 13.2, height = full_h, limitsize = FALSE)
}

plot_pattern_signaling('outgoing', 'SuppFig8B_pattern_outgoing')
plot_pattern_signaling('incoming', 'SuppFig8C_pattern_incoming')
