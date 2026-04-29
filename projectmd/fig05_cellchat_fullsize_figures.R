#!/usr/bin/env Rscript
# fig05_cellchat_fullsize_figures.R
# Full-information, enlarged CellChat figures for the main reproduction line.
# These figures keep the original CellChat content but use larger canvases and
# smaller labels so dense CellChat outputs are readable.

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
    library(NMF)
    library(ggalluvial)
})

out <- output_path('fig5_fullsize')
dir.create(out, showWarnings = FALSE, recursive = TRUE)

cc_path <- output_path('fig5/cellchat_merged_enhanced.rds')
if (!file.exists(cc_path)) {
    cc_path <- output_path('fig5/cellchat_merged.rds')
}
cc <- readRDS(cc_path)

if (is.null(cc@netP$centr)) {
    cc <- netAnalysis_computeCentrality(cc, slot.name = 'netP')
}

if (is.null(cc@netP$pattern$outgoing$pattern) || is.null(cc@netP$pattern$incoming$pattern)) {
    cc <- identifyCommunicationPatterns(cc, pattern = 'outgoing', k = 5)
    cc <- identifyCommunicationPatterns(cc, pattern = 'incoming', k = 5)
}

groups <- levels(cc@idents)
mal_grps <- grep('-like|^Epe$|^NSC$', groups, value = TRUE)
micro_grps <- grep('^Mic|^TC$|^EC$|^Per$', groups, value = TRUE)
group_size <- as.numeric(table(cc@idents))

ct_pal <- c(
    'Ast' = '#E41A1C', 'Ast-like' = '#377EB8', 'EC' = '#4DAF4A',
    'Epe' = '#984EA3', 'Epe-like' = '#FF8C00', 'Mic' = '#F781BF',
    'Neu' = '#B39DDB', 'Neu-like' = '#A65628', 'NSC' = '#56B4E9',
    'NSC-like' = '#253494', 'Oli' = '#1B9E77', 'Oli-like' = '#B2DF8A',
    'OPC' = '#E6C200', 'OPC-like' = '#FB9A99', 'Per' = '#E7298A',
    'RGC-like' = '#8E0152', 'TC' = '#00BFC4'
)
ct_pal <- ct_pal[groups]

message('Groups: ', paste(groups, collapse = ', '))
message('Malignant-like groups: ', paste(mal_grps, collapse = ', '))
message('Microenvironment groups: ', paste(micro_grps, collapse = ', '))

comm_all <- subsetCommunication(cc)
comm_mal_to_micro <- subsetCommunication(cc, sources.use = mal_grps, targets.use = micro_grps)
comm_micro_to_mal <- subsetCommunication(cc, sources.use = micro_grps, targets.use = mal_grps)

lr_label_count <- function(comm) {
    label_col <- if ('interaction_name_2' %in% colnames(comm)) 'interaction_name_2' else 'interaction_name'
    length(unique(comm[[label_col]]))
}

bubble_size <- function(comm, sources, targets) {
    n_x <- length(sources) * length(targets)
    n_y <- lr_label_count(comm)
    c(width = max(24, n_x * 0.48 + 8), height = max(30, n_y * 0.17 + 7))
}

rank_df <- rankNet(cc, mode = 'single', stacked = FALSE, do.stat = FALSE,
                   return.data = TRUE)$signaling.contribution
write.csv(rank_df, file.path(out, 'pathway_importance_fullsize_source.csv'), row.names = FALSE)

# Fig 5B: full all-pair bubble plots, enlarged.
if (length(mal_grps) > 0 && length(micro_grps) > 0) {
    sz <- bubble_size(comm_mal_to_micro, mal_grps, micro_grps)
    pdf(file.path(out, 'Fig5B_bubble_mal_to_micro_fullsize.pdf'), width = sz['width'], height = sz['height'])
    tryCatch(print(netVisual_bubble(
        cc, sources.use = mal_grps, targets.use = micro_grps,
        remove.isolate = TRUE, color.heatmap = 'Spectral',
        font.size = 6.2, font.size.title = 11,
        dot.size.min = 1.2, dot.size.max = 4.0,
        angle.x = 90, line.size = 0.15,
        title.name = 'All ligand-receptor signals: malignant-like -> microenvironment'
    )), error = function(e) message('full bubble mal->micro: ', conditionMessage(e)))
    dev.off()

    sz <- bubble_size(comm_micro_to_mal, micro_grps, mal_grps)
    pdf(file.path(out, 'Fig5B_bubble_micro_to_mal_fullsize.pdf'), width = sz['width'], height = sz['height'])
    tryCatch(print(netVisual_bubble(
        cc, sources.use = micro_grps, targets.use = mal_grps,
        remove.isolate = TRUE, color.heatmap = 'Spectral',
        font.size = 6.2, font.size.title = 11,
        dot.size.min = 1.2, dot.size.max = 4.0,
        angle.x = 90, line.size = 0.15,
        title.name = 'All ligand-receptor signals: microenvironment -> malignant-like'
    )), error = function(e) message('full bubble micro->mal: ', conditionMessage(e)))
    dev.off()
}

# Fig 5B: full crosstalk circles, separated and enlarged.
pdf(file.path(out, 'Fig5B_crosstalk_count_circle_fullsize.pdf'), width = 16, height = 16)
par(xpd = TRUE, mar = c(2, 2, 4, 2))
netVisual_circle(cc@net$count, vertex.weight = group_size, color.use = ct_pal,
                 weight.scale = TRUE, label.edge = FALSE,
                 vertex.label.cex = 1.15, vertex.size.max = 12,
                 edge.width.max = 6, alpha.edge = 0.45,
                 margin = 0.22, title.name = '# interactions')
dev.off()

pdf(file.path(out, 'Fig5B_crosstalk_strength_circle_fullsize.pdf'), width = 16, height = 16)
par(xpd = TRUE, mar = c(2, 2, 4, 2))
netVisual_circle(cc@net$weight, vertex.weight = group_size, color.use = ct_pal,
                 weight.scale = TRUE, label.edge = FALSE,
                 vertex.label.cex = 1.15, vertex.size.max = 12,
                 edge.width.max = 6, alpha.edge = 0.45,
                 margin = 0.22, title.name = 'Interaction strength')
dev.off()

# Supp Fig 8A/D: original summaries, enlarged.
pdf(file.path(out, 'SuppFig8A_interaction_heatmap_fullsize.pdf'), width = 13, height = 11)
print(netVisual_heatmap(cc, measure = 'weight', font.size = 8, font.size.title = 13))
dev.off()

pdf(file.path(out, 'SuppFig8D_signalingRole_scatter_fullsize.pdf'), width = 12, height = 9)
print(netAnalysis_signalingRole_scatter(cc, color.use = ct_pal,
                                        label.size = 4.2, dot.size = c(4, 11),
                                        font.size = 13, font.size.title = 14))
dev.off()

# Supp Fig 8B/C: full NMF river plots, enlarged.
# Native identifyCommunicationPatterns heatmaps put all pathways into a small
# panel and remain unreadable even on large canvases, so the readable heatmaps
# are written by fig05_cellchat_readable_figures.R instead.
tryCatch({
    pdf(file.path(out, 'SuppFig8B_pattern_outgoing_river_fullsize.pdf'), width = 22, height = 28)
    print(netAnalysis_river(cc, pattern = 'outgoing', font.size = 4.2, font.size.title = 16))
    dev.off()

    pdf(file.path(out, 'SuppFig8C_pattern_incoming_river_fullsize.pdf'), width = 22, height = 28)
    print(netAnalysis_river(cc, pattern = 'incoming', font.size = 4.2, font.size.title = 16))
    dev.off()
}, error = function(e) message('full NMF pattern skipped: ', conditionMessage(e)))

# Fig 5C: full pathway circles, enlarged. These keep all edges for the pathway.
focus_pw <- intersect(c('MIF', 'SPP1', 'GALECTIN', 'PTN', 'MK', 'VEGF', 'TGFb', 'CXCL'), cc@netP$pathways)
for (pw in focus_pw) {
    pdf(file.path(out, sprintf('Fig5C_pathway_%s_circle_fullsize.pdf', pw)), width = 13, height = 13)
    par(xpd = TRUE, mar = c(2, 2, 4, 2))
    tryCatch(netVisual_aggregate(cc, signaling = pw, layout = 'circle', color.use = ct_pal,
                                 vertex.label.cex = 1.05, vertex.size.max = 9,
                                 edge.width.max = 6, signaling.name = pw,
                                 title.space = 5),
             error = function(e) message(pw, ': ', conditionMessage(e)))
    dev.off()
}

# A full-size chord panel matching the polished style, but with more room.
chord_pw <- intersect(c('MK', 'MIF', 'PTN', 'SPP1', 'GALECTIN', 'TGFb', 'VEGF'), cc@netP$pathways)
if (length(chord_pw) > 0) {
    nc <- 3
    nr <- ceiling(length(chord_pw) / nc)
    pdf(file.path(out, 'Fig5C_pathway_chord_panel_fullsize.pdf'), width = 7 * nc, height = 7 * nr)
    par(mfrow = c(nr, nc), mar = c(1, 1, 3, 1), xpd = TRUE)
    for (pw in chord_pw) {
        tryCatch(netVisual_aggregate(cc, signaling = pw, layout = 'chord',
                                     color.use = ct_pal, signaling.name = pw,
                                     small.gap = 1, big.gap = 8,
                                     vertex.label.cex = 0.85),
                 error = function(e) { plot.new(); title(paste(pw, 'err')) })
    }
    dev.off()
}

saveRDS(cc, output_path('fig5/cellchat_merged_enhanced.rds'))
message('Done. Full-size figures written to: ', out)