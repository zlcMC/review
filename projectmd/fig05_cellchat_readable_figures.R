#!/usr/bin/env Rscript
# fig05_cellchat_readable_figures.R
# Readability-oriented CellChat figures for Fig 5 / Supp Fig 8.
# Keeps the full outputs in output/fig5 untouched and writes simplified views to
# output/fig5_readable/.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
if (!file.exists('workspace_paths.R') && length(file_arg) > 0) {
    script_file <- normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE)
    setwd(dirname(dirname(script_file)))
}
source('workspace_paths.R')

suppressPackageStartupMessages({
    library(CellChat)
    library(ggplot2)
    library(patchwork)
})

out <- output_path('fig5_readable')
dir.create(out, showWarnings = FALSE, recursive = TRUE)

cc_path <- output_path('fig5/cellchat_merged_enhanced.rds')
if (!file.exists(cc_path)) {
    cc_path <- output_path('fig5/cellchat_merged.rds')
}
cc <- readRDS(cc_path)

groups <- levels(cc@idents)
mal_grps <- grep('-like|^Epe$|^NSC$', groups, value = TRUE)
micro_grps <- grep('^Mic|^TC$|^EC$|^Per$', groups, value = TRUE)

ct_pal <- c(
    'Ast' = '#E41A1C', 'Ast-like' = '#377EB8', 'EC' = '#4DAF4A',
    'Epe' = '#984EA3', 'Epe-like' = '#FF8C00', 'Mic' = '#F781BF',
    'Neu' = '#B39DDB', 'Neu-like' = '#A65628', 'NSC' = '#56B4E9',
    'NSC-like' = '#253494', 'Oli' = '#1B9E77', 'Oli-like' = '#B2DF8A',
    'OPC' = '#E6C200', 'OPC-like' = '#FB9A99', 'Per' = '#E7298A',
    'RGC-like' = '#8E0152', 'TC' = '#00BFC4'
)
ct_pal <- ct_pal[intersect(names(ct_pal), groups)]

message('Groups: ', paste(groups, collapse = ', '))
message('Malignant-like groups: ', paste(mal_grps, collapse = ', '))
message('Microenvironment groups: ', paste(micro_grps, collapse = ', '))

comm_all <- subsetCommunication(cc)

safe_prob_limits <- range(comm_all$prob, na.rm = TRUE)
prob_cols <- c('#315CA8', '#37A7B8', '#86D19E', '#F6D96B', '#E86A45', '#B3133B')

theme_readable <- function(base_size = 8) {
    theme_classic(base_size = base_size) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = base_size),
            axis.title = element_text(size = base_size + 1),
            plot.title = element_text(size = base_size + 3, face = 'bold', hjust = 0),
            legend.title = element_text(size = base_size),
            legend.text = element_text(size = base_size - 1),
            plot.margin = margin(8, 12, 8, 12)
        )
}

plot_direction_bubble <- function(comm, sources, targets, direction_name, title, top_n = 45) {
    source_target_levels <- as.vector(outer(sources, targets, paste, sep = ' -> '))
    comm$source_target <- paste(comm$source, comm$target, sep = ' -> ')
    comm <- comm[comm$source_target %in% source_target_levels, , drop = FALSE]
    if (nrow(comm) == 0) {
        warning('No communication rows for ', direction_name)
        return(invisible(NULL))
    }

    rank_df <- aggregate(prob ~ pathway_name + interaction_name, data = comm, FUN = max)
    rank_df <- rank_df[order(-rank_df$prob, rank_df$pathway_name, rank_df$interaction_name), , drop = FALSE]
    rank_df <- head(rank_df, top_n)
    rank_df$interaction_label <- paste(rank_df$pathway_name, rank_df$interaction_name, sep = ': ')

    comm <- merge(
        comm,
        rank_df[, c('pathway_name', 'interaction_name', 'interaction_label')],
        by = c('pathway_name', 'interaction_name')
    )
    comm$source_target <- factor(comm$source_target, levels = source_target_levels)
    comm$interaction_label <- factor(comm$interaction_label, levels = rev(rank_df$interaction_label))
    comm$p_group <- ifelse(comm$pval < 0.01, 'p < 0.01', '0.01 <= p < 0.05')

    write.csv(
        rank_df,
        file.path(out, paste0('Fig5B_', direction_name, '_top_interactions.csv')),
        row.names = FALSE
    )

    plot_width <- max(12, length(source_target_levels) * 0.22 + 4)
    plot_height <- max(8, length(unique(comm$interaction_label)) * 0.18 + 2)

    p <- ggplot(comm, aes(x = source_target, y = interaction_label)) +
        geom_point(aes(color = prob, size = p_group), alpha = 0.9) +
        scale_color_gradientn(colors = prob_cols, limits = safe_prob_limits, name = 'Communication\nprobability') +
        scale_size_manual(values = c('0.01 <= p < 0.05' = 1.2, 'p < 0.01' = 2.4), name = 'Significance') +
        labs(title = title, x = 'Source -> target cell group', y = 'Pathway: ligand-receptor pair') +
        theme_readable(base_size = 7) +
        theme(panel.grid.major = element_line(color = 'grey88', linewidth = 0.25),
              axis.text.y = element_text(size = 5.7))

    ggsave(file.path(out, paste0('Fig5B_bubble_', direction_name, '_top45.pdf')),
           p, width = plot_width, height = plot_height, limitsize = FALSE)
    invisible(p)
}

comm_mal_to_micro <- subsetCommunication(cc, sources.use = mal_grps, targets.use = micro_grps)
comm_micro_to_mal <- subsetCommunication(cc, sources.use = micro_grps, targets.use = mal_grps)

plot_direction_bubble(
    comm_mal_to_micro, mal_grps, micro_grps,
    'mal_to_micro', 'Top ligand-receptor signals: malignant-like -> microenvironment'
)
plot_direction_bubble(
    comm_micro_to_mal, micro_grps, mal_grps,
    'micro_to_mal', 'Top ligand-receptor signals: microenvironment -> malignant-like'
)

plot_split_bubble <- function(comm, fixed_column, fixed_values, varying_column,
                              direction_name, title_prefix, x_label, top_n = 30) {
    for (fixed_value in fixed_values) {
        filtered_comm <- comm[comm[[fixed_column]] == fixed_value, , drop = FALSE]
        if (nrow(filtered_comm) == 0) {
            next
        }

        rank_table <- aggregate(prob ~ pathway_name + interaction_name, data = filtered_comm, FUN = max)
        rank_table <- rank_table[order(-rank_table$prob, rank_table$pathway_name, rank_table$interaction_name), , drop = FALSE]
        rank_table <- head(rank_table, top_n)
        rank_table$interaction_label <- paste(rank_table$pathway_name, rank_table$interaction_name, sep = ': ')

        filtered_comm <- merge(
            filtered_comm,
            rank_table[, c('pathway_name', 'interaction_name', 'interaction_label')],
            by = c('pathway_name', 'interaction_name')
        )
        filtered_comm$cell_group <- factor(filtered_comm[[varying_column]], levels = groups)
        filtered_comm$interaction_label <- factor(filtered_comm$interaction_label, levels = rev(rank_table$interaction_label))
        filtered_comm$p_group <- ifelse(filtered_comm$pval < 0.01, 'p < 0.01', '0.01 <= p < 0.05')

        output_stem <- paste0('Fig5B_bubble_', direction_name, '_', fixed_value, '_top30')
        write.csv(rank_table, file.path(out, paste0(output_stem, '.csv')), row.names = FALSE)

        split_plot <- ggplot(filtered_comm, aes(x = cell_group, y = interaction_label)) +
            geom_point(aes(color = prob, size = p_group), alpha = 0.9) +
            scale_color_gradientn(colors = prob_cols, limits = safe_prob_limits, name = 'Communication\nprobability') +
            scale_size_manual(values = c('0.01 <= p < 0.05' = 1.3, 'p < 0.01' = 2.8), name = 'Significance') +
            labs(title = paste0(title_prefix, fixed_value), x = x_label, y = 'Pathway: ligand-receptor pair') +
            theme_readable(base_size = 8) +
            theme(panel.grid.major = element_line(color = 'grey88', linewidth = 0.25),
                  axis.text.y = element_text(size = 6.3),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

        ggsave(file.path(out, paste0(output_stem, '.pdf')), split_plot,
               width = 8.2, height = max(7, length(unique(filtered_comm$interaction_label)) * 0.19 + 2),
               limitsize = FALSE)
    }
}

plot_split_bubble(
    comm_mal_to_micro, 'target', micro_grps, 'source',
    'mal_to_micro_target', 'Malignant-like -> ', 'Source malignant-like cell group'
)
plot_split_bubble(
    comm_micro_to_mal, 'source', micro_grps, 'target',
    'micro_to_mal_source', 'Microenvironment source: ', 'Target malignant-like cell group'
)

# Overall crosstalk: use heatmaps and top-edge bars instead of dense all-edge circles.
plot_matrix_heatmap <- function(mat, title, file, width = 8.5, height = 7.5) {
    df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(df) <- c('source', 'target', 'strength')
    df$source <- factor(df$source, levels = rev(rownames(mat)))
    df$target <- factor(df$target, levels = colnames(mat))
    p <- ggplot(df, aes(x = target, y = source, fill = strength)) +
        geom_tile(color = 'white', linewidth = 0.25) +
        scale_fill_gradientn(colors = c('#FFF7EC', '#FDD49E', '#FC8D59', '#D7301F', '#7F0000'), name = 'Strength') +
        labs(title = title, x = 'Receiver / target', y = 'Sender / source') +
        coord_fixed() +
        theme_readable(base_size = 7) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.border = element_rect(color = 'grey50', fill = NA, linewidth = 0.3))
    ggsave(file.path(out, file), p, width = width, height = height, limitsize = FALSE)
    invisible(p)
}

plot_top_edges <- function(mat, title, file, top_n = 40) {
    df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(df) <- c('source', 'target', 'strength')
    df <- df[df$strength > 0, , drop = FALSE]
    df <- df[order(-df$strength), , drop = FALSE]
    top_df <- head(df, top_n)
    top_df$edge <- paste(top_df$source, top_df$target, sep = ' -> ')
    top_df$edge <- factor(top_df$edge, levels = rev(top_df$edge))
    write.csv(top_df, file.path(out, sub('pdf$', 'csv', file)), row.names = FALSE)
    p <- ggplot(top_df, aes(x = strength, y = edge, fill = source)) +
        geom_col(width = 0.78, show.legend = FALSE) +
        scale_fill_manual(values = ct_pal, na.value = '#999999') +
        labs(title = title, x = 'Aggregated interaction strength', y = 'Source -> target') +
        theme_classic(base_size = 8) +
        theme(plot.title = element_text(size = 11, face = 'bold'),
              axis.text.y = element_text(size = 6.5),
              plot.margin = margin(8, 12, 8, 12))
    ggsave(file.path(out, file), p, width = 8.3, height = 9.2, limitsize = FALSE)
    invisible(p)
}

plot_matrix_heatmap(cc@net$weight, 'Overall CellChat interaction strength', 'Fig5B_crosstalk_strength_heatmap_readable.pdf')
plot_matrix_heatmap(cc@net$count, 'Overall CellChat interaction count', 'Fig5B_crosstalk_count_heatmap_readable.pdf')
plot_top_edges(cc@net$weight, 'Top source-target pairs by interaction strength', 'Fig5B_crosstalk_top40_edges.pdf')

# Pathway-level heatmaps are more legible than all-edge circle/chord plots for dense pathways.
plot_pathway_heatmap <- function(pathway, individual = TRUE) {
    mat <- cc@netP$prob[, , pathway]
    df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(df) <- c('source', 'target', 'prob')
    df$source <- factor(df$source, levels = rev(rownames(mat)))
    df$target <- factor(df$target, levels = colnames(mat))
    p <- ggplot(df, aes(x = target, y = source, fill = prob)) +
        geom_tile(color = 'white', linewidth = 0.25) +
        scale_fill_gradientn(colors = prob_cols, name = 'Prob.') +
        labs(title = paste0(pathway, ' pathway'), x = 'Receiver / target', y = 'Sender / source') +
        coord_fixed() +
        theme_readable(base_size = if (individual) 7 else 5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = if (individual) 11 else 7, face = 'bold'),
              legend.position = if (individual) 'right' else 'none')
    if (individual) {
        ggsave(file.path(out, paste0('Fig5C_pathway_heatmap_', pathway, '.pdf')),
               p, width = 8.5, height = 7.5, limitsize = FALSE)
    }
    p
}

focus_pw <- intersect(c('MK', 'PTN', 'MIF', 'SPP1', 'GALECTIN', 'TGFb', 'VEGF'), cc@netP$pathways)
pathway_plots <- lapply(focus_pw, function(pathway) plot_pathway_heatmap(pathway, individual = TRUE))
names(pathway_plots) <- focus_pw
panel <- wrap_plots(lapply(focus_pw, function(pathway) plot_pathway_heatmap(pathway, individual = FALSE)), ncol = 4)
ggsave(file.path(out, 'Fig5C_pathway_heatmap_panel_readable.pdf'), panel, width = 15, height = 8.5, limitsize = FALSE)

# NMF pattern summaries: top signaling per pattern and cell-group contribution heatmaps.
plot_pattern_signaling <- function(direction, file_prefix, top_n = 10) {
    pat <- cc@netP$pattern[[direction]]$pattern
    if (is.null(pat) || is.null(pat$signaling) || is.null(pat$cell)) {
        warning('No pattern data for ', direction)
        return(invisible(NULL))
    }
    sig <- pat$signaling
    sig$Pattern <- as.character(sig$Pattern)
    sig$Signaling <- as.character(sig$Signaling)
    idx <- unlist(lapply(split(seq_len(nrow(sig)), sig$Pattern), function(i) {
        i <- i[order(-sig$Contribution[i])]
        head(i, top_n)
    }), use.names = FALSE)
    sig_top <- sig[idx, , drop = FALSE]
    sig_top <- sig_top[order(sig_top$Pattern, -sig_top$Contribution), , drop = FALSE]
    sig_top$Signaling <- factor(sig_top$Signaling, levels = rev(unique(sig_top$Signaling)))
    write.csv(sig_top, file.path(out, paste0(file_prefix, '_top_signaling.csv')), row.names = FALSE)

    p_sig <- ggplot(sig_top, aes(x = Pattern, y = Signaling)) +
        geom_point(aes(size = Contribution, color = Contribution), alpha = 0.92) +
        scale_color_gradientn(colors = prob_cols, name = 'Contribution') +
        scale_size(range = c(1.2, 4.2), name = 'Contribution') +
        labs(title = paste0(direction, ' patterns: top signaling pathways'), x = NULL, y = 'Signaling pathway') +
        theme_classic(base_size = 8) +
        theme(plot.title = element_text(size = 11, face = 'bold'),
              axis.text.y = element_text(size = 6.5),
              panel.grid.major.y = element_line(color = 'grey88', linewidth = 0.25))
    ggsave(file.path(out, paste0(file_prefix, '_top_signaling_dotplot.pdf')),
           p_sig, width = 7.5, height = max(6, length(unique(sig_top$Signaling)) * 0.17 + 2), limitsize = FALSE)

    sig_full <- sig
    dominant_idx <- unlist(lapply(split(seq_len(nrow(sig_full)), sig_full$Signaling), function(i) {
        i[which.max(sig_full$Contribution[i])]
    }), use.names = FALSE)
    dominant <- sig_full[dominant_idx, , drop = FALSE]
    dominant <- dominant[order(dominant$Pattern, -dominant$Contribution, dominant$Signaling), , drop = FALSE]
    sig_full$Pattern <- factor(sig_full$Pattern, levels = sort(unique(sig_full$Pattern)))
    sig_full$Signaling <- factor(sig_full$Signaling, levels = rev(dominant$Signaling))
    write.csv(sig_full, file.path(out, paste0(file_prefix, '_all_signaling.csv')), row.names = FALSE)

    p_sig_full <- ggplot(sig_full, aes(x = Pattern, y = Signaling, fill = Contribution)) +
        geom_tile(color = 'white', linewidth = 0.18) +
        scale_fill_gradientn(colors = prob_cols, name = 'Contribution') +
        labs(title = paste0(direction, ' patterns: all signaling pathways'), x = NULL, y = 'Signaling pathway') +
        theme_classic(base_size = 8) +
        theme(plot.title = element_text(size = 11, face = 'bold'),
              axis.text.x = element_text(angle = 0),
              axis.text.y = element_text(size = 5.2),
              panel.border = element_rect(color = 'grey50', fill = NA, linewidth = 0.3))
    full_height <- max(11, length(unique(sig_full$Signaling)) * 0.13 + 2)
    ggsave(file.path(out, paste0(file_prefix, '_all_signaling_heatmap.pdf')),
           p_sig_full, width = 6.6, height = full_height, limitsize = FALSE)

    cell <- pat$cell
    cell$Pattern <- as.character(cell$Pattern)
    cell$CellGroup <- factor(as.character(cell$CellGroup), levels = rev(groups))
    p_cell <- ggplot(cell, aes(x = Pattern, y = CellGroup, fill = Contribution)) +
        geom_tile(color = 'white', linewidth = 0.25) +
        scale_fill_gradientn(colors = prob_cols, name = 'Contribution') +
        labs(title = paste0(direction, ' patterns: cell-group contributions'), x = NULL, y = NULL) +
        theme_classic(base_size = 8) +
        theme(plot.title = element_text(size = 11, face = 'bold'),
              axis.text.x = element_text(angle = 0),
              panel.border = element_rect(color = 'grey50', fill = NA, linewidth = 0.3))
    ggsave(file.path(out, paste0(file_prefix, '_cell_contribution_heatmap.pdf')),
           p_cell, width = 5.8, height = 5.8, limitsize = FALSE)

        overview <- p_cell + p_sig_full + plot_layout(widths = c(0.95, 1.35))
        ggsave(file.path(out, paste0(file_prefix, '_readable_overview.pdf')),
            overview, width = 13.2, height = full_height, limitsize = FALSE)
}

plot_pattern_signaling('outgoing', 'SuppFig8B_pattern_outgoing')
plot_pattern_signaling('incoming', 'SuppFig8C_pattern_incoming')

message('Done. Readable figures written to: ', out)