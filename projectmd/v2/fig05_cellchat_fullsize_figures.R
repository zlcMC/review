#!/usr/bin/env Rscript
# v2 of fig05_cellchat_fullsize_figures.R
# Full-information enlarged CellChat figures.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(CellChat); library(ggplot2); library(patchwork)
  library(NMF); library(ggalluvial)
})

out <- v2_dir('fig5_fullsize')
cc  <- cellchat_load()

if (is.null(cc@netP$centr)) {
  cc <- netAnalysis_computeCentrality(cc, slot.name = 'netP')
}
if (is.null(cc@netP$pattern$outgoing$pattern) ||
    is.null(cc@netP$pattern$incoming$pattern)) {
  cc <- identifyCommunicationPatterns(cc, pattern = 'outgoing', k = 5)
  cc <- identifyCommunicationPatterns(cc, pattern = 'incoming', k = 5)
}

g       <- cellchat_groups(cc)
ct_pal  <- cellchat_palette(g$all, restrict = FALSE)[g$all]
group_size <- as.numeric(table(cc@idents))

message('Groups: ', paste(g$all, collapse = ', '))
message('Malignant-like: ', paste(g$mal, collapse = ', '))
message('Microenvironment: ', paste(g$micro, collapse = ', '))

comm_mal_to_micro <- subsetCommunication(cc, sources.use = g$mal, targets.use = g$micro)
comm_micro_to_mal <- subsetCommunication(cc, sources.use = g$micro, targets.use = g$mal)

lr_label_count <- function(comm) {
  col <- if ('interaction_name_2' %in% colnames(comm)) 'interaction_name_2' else 'interaction_name'
  length(unique(comm[[col]]))
}
bubble_size <- function(comm, sources, targets) {
  c(width  = max(24, length(sources) * length(targets) * 0.48 + 8),
    height = max(30, lr_label_count(comm) * 0.17 + 7))
}

rank_df <- rankNet(cc, mode = 'single', stacked = FALSE, do.stat = FALSE,
                   return.data = TRUE)$signaling.contribution
write.csv(rank_df, file.path(out, 'pathway_importance_fullsize_source.csv'), row.names = FALSE)

draw_full_bubble <- function(sources, targets, comm, file, title) {
  if (!length(sources) || !length(targets) || !nrow(comm)) return(invisible())
  sz <- bubble_size(comm, sources, targets)
  pdf(file.path(out, file), width = sz['width'], height = sz['height'])
  tryCatch(print(netVisual_bubble(
    cc, sources.use = sources, targets.use = targets,
    remove.isolate = TRUE, color.heatmap = 'Spectral',
    font.size = 6.2, font.size.title = 11,
    dot.size.min = 1.2, dot.size.max = 4.0,
    angle.x = 90, line.size = 0.15, title.name = title
  )), error = function(e) message(file, ': ', conditionMessage(e)))
  dev.off()
}
draw_full_bubble(g$mal, g$micro, comm_mal_to_micro,
                 'Fig5B_bubble_mal_to_micro_fullsize.pdf',
                 'All ligand-receptor signals: malignant-like -> microenvironment')
draw_full_bubble(g$micro, g$mal, comm_micro_to_mal,
                 'Fig5B_bubble_micro_to_mal_fullsize.pdf',
                 'All ligand-receptor signals: microenvironment -> malignant-like')

draw_circle <- function(mat, file, title) {
  pdf(file.path(out, file), width = 16, height = 16)
  par(xpd = TRUE, mar = c(2, 2, 4, 2))
  netVisual_circle(mat, vertex.weight = group_size, color.use = ct_pal,
                   weight.scale = TRUE, label.edge = FALSE,
                   vertex.label.cex = 1.15, vertex.size.max = 12,
                   edge.width.max = 6, alpha.edge = 0.45,
                   margin = 0.22, title.name = title)
  dev.off()
}
draw_circle(cc@net$count,  'Fig5B_crosstalk_count_circle_fullsize.pdf',    '# interactions')
draw_circle(cc@net$weight, 'Fig5B_crosstalk_strength_circle_fullsize.pdf', 'Interaction strength')

pdf(file.path(out, 'SuppFig8A_interaction_heatmap_fullsize.pdf'), width = 13, height = 11)
print(netVisual_heatmap(cc, measure = 'weight', font.size = 8, font.size.title = 13))
dev.off()

pdf(file.path(out, 'SuppFig8D_signalingRole_scatter_fullsize.pdf'), width = 12, height = 9)
print(netAnalysis_signalingRole_scatter(cc, color.use = ct_pal,
      label.size = 4.2, dot.size = c(4, 11), font.size = 13, font.size.title = 14))
dev.off()

tryCatch({
  pdf(file.path(out, 'SuppFig8B_pattern_outgoing_river_fullsize.pdf'), width = 22, height = 28)
  print(netAnalysis_river(cc, pattern = 'outgoing', font.size = 4.2, font.size.title = 16))
  dev.off()
  pdf(file.path(out, 'SuppFig8C_pattern_incoming_river_fullsize.pdf'), width = 22, height = 28)
  print(netAnalysis_river(cc, pattern = 'incoming', font.size = 4.2, font.size.title = 16))
  dev.off()
}, error = function(e) message('full NMF pattern skipped: ', conditionMessage(e)))

focus_pw <- intersect(c('MIF','SPP1','GALECTIN','PTN','MK','VEGF','TGFb','CXCL'), cc@netP$pathways)
for (pw in focus_pw) {
  pdf(file.path(out, sprintf('Fig5C_pathway_%s_circle_fullsize.pdf', pw)), width = 13, height = 13)
  par(xpd = TRUE, mar = c(2, 2, 4, 2))
  tryCatch(netVisual_aggregate(cc, signaling = pw, layout = 'circle', color.use = ct_pal,
                               vertex.label.cex = 1.05, vertex.size.max = 9,
                               edge.width.max = 6, signaling.name = pw, title.space = 5),
           error = function(e) message(pw, ': ', conditionMessage(e)))
  dev.off()
}

chord_pw <- intersect(c('MK','MIF','PTN','SPP1','GALECTIN','TGFb','VEGF'), cc@netP$pathways)
if (length(chord_pw)) {
  nc <- 3; nr <- ceiling(length(chord_pw) / nc)
  pdf(file.path(out, 'Fig5C_pathway_chord_panel_fullsize.pdf'), width = 7 * nc, height = 7 * nr)
  par(mfrow = c(nr, nc), mar = c(1, 1, 3, 1), xpd = TRUE)
  for (pw in chord_pw) {
    tryCatch(netVisual_aggregate(cc, signaling = pw, layout = 'chord',
                                 color.use = ct_pal, signaling.name = pw,
                                 small.gap = 1, big.gap = 8, vertex.label.cex = 0.85),
             error = function(e) { plot.new(); title(paste(pw, 'err')) })
  }
  dev.off()
}

saveRDS(cc, v2_file('fig5', 'cellchat_merged_enhanced.rds'))
