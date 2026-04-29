#!/usr/bin/env bash
# Run all v2 R scripts in fig order, log per-script status.
set -u
cd "$(dirname "$0")/../.."
LOG_DIR="output/v2/_runtime_logs"
mkdir -p "$LOG_DIR"
SUMMARY="$LOG_DIR/summary.tsv"
: > "$SUMMARY"

SCRIPTS=(
  cohort_prepare_pajtler_template.R
  fig01_bulk_relapse_deg.R
  fig01_supp_subclone_deg_go.R
  fig01_tsne_cnv_cytotrace_seurat.R
  fig02_slingshot_basic.R
  fig02c_slingshot_final_panel.R
  fig03_subclone_undiff_diff_score.R
  fig04_bulk_relapse_proxies.R
  fig04_gillen_microarray_proxy.R
  fig04_gillen_scrna_deg_go.R
  fig04_gojo_scrna_deg_go.R
  fig04_gillen_survival_km.R
  fig04_survival_template.R
  fig05_cellchat_prepare.R
  fig05_cellchat_downstream.R
  fig05_cellchat_fullsize_figures.R
  fig05_cellchat_readable_figures.R
  final_paper_style_export.R
)

TIMEOUT="${TIMEOUT:-600}"
echo "timeout per script: ${TIMEOUT}s"
printf 'status\tdur(s)\tscript\n' | tee "$SUMMARY"

for s in "${SCRIPTS[@]}"; do
  log="$LOG_DIR/${s%.R}.log"
  start=$(date +%s)
  timeout "$TIMEOUT" conda run -n epn2_r --no-capture-output \
    Rscript "projectmd/v2/$s" > "$log" 2>&1
  rc=$?
  dur=$(( $(date +%s) - start ))
  case $rc in
    0)   st='PASS' ;;
    124) st='TIMEOUT' ;;
    *)   st="FAIL($rc)" ;;
  esac
  printf '%s\t%s\t%s\n' "$st" "$dur" "$s" | tee -a "$SUMMARY"
done

echo
echo '== summary =='
column -t "$SUMMARY"
