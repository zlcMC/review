#!/usr/bin/env bash
# fig06_scenic_ctx_aucell_pipeline.sh — 在 GRN 完成后跑 ctx + aucell + 下游
set -e
cd "$(dirname "$0")/.."

LOOM=output/fig6/scenic_input/integrated.loom
ADJ=output/fig6/adj.tsv
CIS=output/fig6/cistarget
REG=output/fig6/regulons.csv
AUC=output/fig6/auc_mtx.loom

[[ -s "$ADJ" ]] || { echo "ERROR: $ADJ empty/missing — GRN not done"; exit 1; }
echo "✓ adj.tsv: $(wc -l < $ADJ) lines"

# ctx: prune to direct targets via cisTarget motif DB
FEATHER=$(ls $CIS/*.feather | head -1)
MOTIFS=$(ls $CIS/*motifs*.tbl 2>/dev/null | head -1)
echo "Feather: $FEATHER"
echo "Motifs: $MOTIFS"

echo '=== ctx ==='
conda run -n epn2 pyscenic ctx "$ADJ" "$FEATHER" \
    --annotations_fname "$MOTIFS" \
    --expression_mtx_fname "$LOOM" \
    --output "$REG" --mask_dropouts --num_workers 4

echo '=== aucell ==='
conda run -n epn2 pyscenic aucell "$LOOM" "$REG" \
    --output "$AUC" --num_workers 4

echo '=== downstream Fig 6 ==='
conda run -n epn2 python projectmd/fig06_pyscenic_downstream.py

echo '✅ pySCENIC pipeline + Fig 6 done'
