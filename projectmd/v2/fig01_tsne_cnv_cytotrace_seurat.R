#!/usr/bin/env Rscript
# v2 of fig01_tsne_cnv_cytotrace_seurat.R
# GTE001/002/009/012: counts → QC → HVG → CC regression → PCA/tSNE → cluster.

this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
})
options(Seurat.object.assay.version = 'v3')

SAMPLES <- c('GTE001', 'GTE002', 'GTE009', 'GTE012')
GSM <- c(GTE001 = 'GSM5710226', GTE002 = 'GSM5710227',
         GTE009 = 'GSM5710228', GTE012 = 'GSM5710229')

process_one <- function(s) {
  cat('\n====', s, '====\n')
  cache <- v2_file('fig1', paste0(s, '_seurat.rds'))
  fallback <- dep_file('fig1', paste0(s, '_seurat.rds'))
  if (file.exists(cache))   { cat('  (v2 cache)\n');  return(readRDS(cache)) }
  if (file.exists(fallback) && fallback != cache) {
    cat('  (legacy cache)\n'); return(readRDS(fallback))
  }
  t0 <- Sys.time()
  df   <- read.csv(raw_data_path(paste0(GSM[s], '_', s, '_counts.csv.gz')),
                   row.names = 1, check.names = FALSE)
  meta <- read.csv(raw_data_path(paste0(GSM[s], '_', s, '_metadata.csv.gz')),
                   row.names = 1, check.names = FALSE)
  cat(sprintf('  read: %d × %d  (%.1fs)\n', nrow(df), ncol(df),
              as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  obj <- CreateSeuratObject(counts = df, meta.data = meta, min.cells = 10)
  rm(df); invisible(gc())

  if (!'percent.mt' %in% colnames(obj@meta.data)) {
    obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = '^MT-')
  }
  obj <- subset(obj, subset = nFeature_RNA >= 1500 & percent.mt < 12)
  if ('DF_hi.lo' %in% colnames(obj@meta.data)) {
    obj <- subset(obj, subset = DF_hi.lo == 'Singlet')
  }
  cat(sprintf('  QC+Singlet: %d cells\n', ncol(obj)))

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 5000, verbose = FALSE)
  if (!all(c('S.Score', 'G2M.Score') %in% colnames(obj@meta.data))) {
    obj <- CellCycleScoring(obj,
            s.features  = cc.genes.updated.2019$s.genes,
            g2m.features = cc.genes.updated.2019$g2m.genes,
            set.ident = FALSE)
  }
  obj <- ScaleData(obj, vars.to.regress = c('S.Score', 'G2M.Score'), verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- RunTSNE(obj, dims = 1:50, seed.use = 42)
  obj <- FindNeighbors(obj, dims = 1:50, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  cat(sprintf('  done in %.1fs\n', as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  saveRDS(obj, cache)
  cat('  cached →', cache, '\n')
  obj
}

plots_for_sample <- function(obj, s) {
  spacer <- patchwork::plot_spacer()
  p_bc  <- DimPlot(obj, reduction = 'tsne', group.by = 'Brief_cluster',
                   label = TRUE, pt.size = 0.3) + ggtitle(paste(s, 'Brief_cluster'))
  p_cnv <- if ('CNV_level' %in% colnames(obj@meta.data))
    FeaturePlot(obj, reduction = 'tsne', features = 'CNV_level',
                cols = c('grey90', 'firebrick')) +
      ggtitle(paste(s, 'CNV level (Fig 1A)'))
    else spacer
  p_cyt <- if ('CytoTRACE' %in% colnames(obj@meta.data))
    FeaturePlot(obj, reduction = 'tsne', features = 'CytoTRACE',
                cols = c('blue', 'white', 'red')) +
      ggtitle(paste(s, 'CytoTRACE (Fig 1B)'))
    else spacer
  p_sub <- if ('CNV_Cluster' %in% colnames(obj@meta.data))
    DimPlot(obj, reduction = 'tsne', group.by = 'CNV_Cluster',
            label = TRUE, pt.size = 0.3) +
      ggtitle(paste(s, 'Subclone (Fig 1C/D)'))
    else spacer
  out <- v2_file('fig1', paste0(s, '_tsne_panels.pdf'))
  ggsave(out, (p_bc | p_cnv) / (p_cyt | p_sub), width = 14, height = 12)
  cat('  fig →', out, '\n')
}

args <- commandArgs(trailingOnly = TRUE)
samples <- if (length(args)) args else SAMPLES
for (s in samples) {
  obj <- process_one(s)
  plots_for_sample(obj, s)
}
cat('\n✓ all done\n')
