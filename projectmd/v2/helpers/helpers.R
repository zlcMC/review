# projectmd/v2/helpers/helpers.R
# Shared R helpers for the v2 refactored scripts.
#
# Source from any v2 R script via the one-liner:
#   source(file.path(dirname(sys.frame(1)$ofile), 'helpers/helpers.R'))
# or, more robustly, via the wrapper init_v2() defined at the bottom which
# both sets up the workspace and sources this file.

# ---------------------------------------------------------------------------
# 1) Workspace bootstrap
# ---------------------------------------------------------------------------

init_workspace <- function() {
  if (file.exists('workspace_paths.R')) {
    source('workspace_paths.R', chdir = FALSE)
    return(invisible(TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep('^--file=', args, value = TRUE)
  if (!length(fa)) {
    stop('Run this script from the workspace root, or via `Rscript projectmd/v2/...R`.')
  }
  sf <- normalizePath(sub('^--file=', '', fa[1]), winslash = '/', mustWork = TRUE)
  repo <- dirname(sf)
  while (!file.exists(file.path(repo, 'workspace_paths.R')) && dirname(repo) != repo) {
    repo <- dirname(repo)
  }
  if (!file.exists(file.path(repo, 'workspace_paths.R'))) {
    stop('Could not locate workspace_paths.R from ', sf)
  }
  setwd(repo)
  source('workspace_paths.R', chdir = FALSE)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# 2) v2 output / dependency path helpers
# ---------------------------------------------------------------------------

v2_file <- function(...) {
  p <- file.path(output_dir, 'v2', ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

v2_dir <- function(...) {
  p <- file.path(output_dir, 'v2', ...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# Read a dependency: prefer output/v2/<...>, fall back to output/<...>.
dep_file <- function(...) {
  v2 <- file.path(output_dir, 'v2', ...)
  if (file.exists(v2)) return(v2)
  file.path(output_dir, ...)
}

# ---------------------------------------------------------------------------
# 3) Generic GEO series-matrix loader
#    Replaces the copy-pasted load_gse64415 / load_gse125861 blocks across
#    fig01_bulk_relapse_deg.R, fig04_bulk_relapse_proxies.R,
#    fig04_gillen_microarray_proxy.R and fig04_gillen_survival_km.R.
# ---------------------------------------------------------------------------

load_gse_expression <- function(matrix_path,
                                extra_columns = list(),
                                sample_filter_key = NULL,
                                sample_filter_values = NULL,
                                annotate_probes = TRUE,
                                require_log2 = TRUE,
                                numeric_columns = c('age')) {
  suppressPackageStartupMessages({
    library(GEOquery)
    library(limma)
  })
  gse  <- getGEO(filename = matrix_path, getGPL = FALSE)
  expr <- exprs(gse)
  pd   <- pData(gse)
  chars <- pd[, grep('characteristics', colnames(pd)), drop = FALSE]
  extract_field <- function(key) {
    apply(chars, 1, function(row) {
      hit <- grep(key, row, value = TRUE, ignore.case = TRUE)
      if (length(hit)) sub('.*: *', '', hit[1]) else NA_character_
    })
  }
  meta <- data.frame(
    sample = rownames(pd),
    title  = pd$title,
    stringsAsFactors = FALSE
  )
  for (nm in names(extra_columns)) {
    val <- extract_field(extra_columns[[nm]])
    if (nm %in% numeric_columns) val <- suppressWarnings(as.numeric(val))
    meta[[nm]] <- val
  }
  if (!is.null(sample_filter_key)) {
    keep <- !is.na(meta[[sample_filter_key]]) &
            meta[[sample_filter_key]] %in% sample_filter_values
    expr <- expr[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
  }
  if (require_log2 && length(expr) && max(expr, na.rm = TRUE) > 50) {
    expr <- log2(expr + 1)
  }
  expr <- expr[rowSums(is.na(expr)) == 0, , drop = FALSE]
  if (annotate_probes && nrow(expr) > 0 && grepl('_at$', rownames(expr)[1])) {
    suppressPackageStartupMessages(library(hgu133plus2.db))
    map <- suppressMessages(AnnotationDbi::select(
      hgu133plus2.db,
      keys     = rownames(expr),
      columns  = 'SYMBOL',
      keytype  = 'PROBEID'
    ))
    map  <- map[!is.na(map$SYMBOL) & !duplicated(map$PROBEID), ]
    expr <- expr[map$PROBEID, , drop = FALSE]
    expr <- limma::avereps(expr, ID = map$SYMBOL)
  }
  list(expr = expr, meta = meta)
}

# ---------------------------------------------------------------------------
# 4) Trajectory signature scoring
#    Replaces score_bulk_trajectory / score_with_signatures across four scripts.
# ---------------------------------------------------------------------------

load_trajectory_signatures <- function(
    rds_path = dep_file('fig3', 'trajectory_signatures.rds')) {
  if (!file.exists(rds_path)) {
    stop('Missing trajectory signatures: ', rds_path)
  }
  sigs <- readRDS(rds_path)
  list(
    undiff = if (!is.null(sigs$undiff_genes)) sigs$undiff_genes else sigs$undiff,
    diff   = if (!is.null(sigs$diff_genes))   sigs$diff_genes   else sigs$diff
  )
}

score_trajectory <- function(expr, meta,
                             sample_col = 'sample',
                             min_overlap = 5,
                             rds_path = dep_file('fig3', 'trajectory_signatures.rds')) {
  sigs   <- load_trajectory_signatures(rds_path)
  undiff <- intersect(sigs$undiff, rownames(expr))
  diff   <- intersect(sigs$diff,   rownames(expr))
  if (length(undiff) < min_overlap || length(diff) < min_overlap) {
    stop(sprintf('Too few signature genes matched (undiff=%d, diff=%d).',
                 length(undiff), length(diff)))
  }
  score <- colMeans(expr[undiff, , drop = FALSE], na.rm = TRUE) -
           colMeans(expr[diff,   , drop = FALSE], na.rm = TRUE)
  meta$trajectory_score <- as.numeric(score[meta[[sample_col]]])
  thr <- median(meta$trajectory_score, na.rm = TRUE)
  meta$trajectory_group <- factor(
    ifelse(meta$trajectory_score >= thr, 'High', 'Low'),
    levels = c('Low', 'High')
  )
  attr(meta, 'matched_undiff') <- length(undiff)
  attr(meta, 'matched_diff')   <- length(diff)
  meta
}

# ---------------------------------------------------------------------------
# 5) GO enrichment + dotplot helper
#    Replaces the duplicated enrichGO / dotplot / ggsave block across
#    fig04_gillen_scrna_deg_go.R and fig04_gojo_scrna_deg_go.R.
#    Returns the enrichResult object invisibly.
# ---------------------------------------------------------------------------

enrich_go_dotplot <- function(genes, title, out_pdf, out_png = NULL,
                              ont = 'BP', pAdjust = 'BH',
                              pcutoff = 0.1, qcutoff = 0.2,
                              show_n = 15, width = 7.4, height = 5.4,
                              csv_path = NULL, genes_txt = NULL,
                              csv_writer = c('readr', 'base')) {
  csv_writer <- match.arg(csv_writer)
  suppressPackageStartupMessages({
    library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
  })
  genes <- unique(genes)
  if (!is.null(genes_txt)) writeLines(genes, genes_txt)
  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                  ont = ont, pAdjustMethod = pAdjust,
                  pvalueCutoff = pcutoff, qvalueCutoff = qcutoff,
                  readable = TRUE)
  ego_df <- as.data.frame(ego)
  if (!is.null(csv_path)) {
    if (csv_writer == 'base') write.csv(ego_df, csv_path, row.names = FALSE)
    else readr::write_csv(ego_df, csv_path)
  }
  if (nrow(ego_df) > 0) {
    p <- dotplot(ego, showCategory = min(show_n, nrow(ego_df))) +
         ggtitle(title) + theme_bw(base_size = 9)
    ggsave(out_pdf, p, width = width, height = height)
    if (!is.null(out_png)) ggsave(out_png, p, width = width, height = height, dpi = 220)
  } else {
    cat('No enriched GO terms passed cutoffs (', title, ')\n')
  }
  cat(sprintf('GO input genes: %d; enriched terms: %d\n', length(genes), nrow(ego_df)))
  invisible(ego)
}

# ---------------------------------------------------------------------------
# 6) CellChat helpers
#    Shared by fig05_cellchat_*.R: object loading, cell-type palette and
#    canonical malignant / microenvironment group splits.
# ---------------------------------------------------------------------------

CELLCHAT_PALETTE <- c(
  'Ast' = '#E41A1C', 'Ast-like' = '#377EB8', 'EC' = '#4DAF4A',
  'Epe' = '#984EA3', 'Epe-like' = '#FF8C00', 'Mic' = '#F781BF',
  'Neu' = '#B39DDB', 'Neu-like' = '#A65628', 'NSC' = '#56B4E9',
  'NSC-like' = '#253494', 'Oli' = '#1B9E77', 'Oli-like' = '#B2DF8A',
  'OPC' = '#E6C200', 'OPC-like' = '#FB9A99', 'Per' = '#E7298A',
  'RGC-like' = '#8E0152', 'TC' = '#00BFC4'
)

cellchat_palette <- function(groups, restrict = TRUE) {
  pal <- CELLCHAT_PALETTE
  if (restrict) pal <- pal[intersect(names(pal), groups)]
  pal
}

cellchat_load <- function(enhanced_first = TRUE) {
  paths <- c('fig5/cellchat_merged_enhanced.rds', 'fig5/cellchat_merged.rds')
  if (!enhanced_first) paths <- rev(paths)
  for (p in paths) {
    full <- dep_file(p)
    if (file.exists(full)) return(readRDS(full))
  }
  stop('No cellchat object found; run fig05_cellchat_prepare.R first.')
}

cellchat_groups <- function(cc) {
  groups <- levels(cc@idents)
  list(
    all   = groups,
    mal   = grep('-like|^Epe$|^NSC$', groups, value = TRUE),
    micro = grep('^Mic|^TC$|^EC$|^Per$', groups, value = TRUE)
  )
}

# ---------------------------------------------------------------------------
# locate this helpers.R and then call init_workspace().
#   this_dir <- (function() {
#     a <- commandArgs(trailingOnly = FALSE)
#     dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
#   })()
#   source(file.path(this_dir, 'helpers', 'helpers.R'))
#   init_workspace()
# ---------------------------------------------------------------------------
