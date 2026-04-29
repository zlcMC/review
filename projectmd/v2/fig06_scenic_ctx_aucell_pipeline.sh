#!/usr/bin/env bash
# v2 of projectmd/fig06_scenic_ctx_aucell_pipeline.sh — ctx + AUCell + 下游
# 写入 output/v2/fig6/，cisTarget 数据库优先 v2，缺失回退 v1。
set -e
cd "$(dirname "$0")/../.."

V2=output/v2/fig6
V1=output/fig6
LOOM=$V2/scenic_input/integrated.loom
ADJ=$V2/adj.tsv
[[ -s "$ADJ" ]] || ADJ=$V1/adj.tsv
[[ -s "$ADJ" ]] || { echo "ERROR: adj.tsv missing in $V2 and $V1"; exit 1; }
[[ -s "$LOOM" ]] || LOOM=$V1/scenic_input/integrated.loom
[[ -s "$LOOM" ]] || { echo "ERROR: integrated.loom missing"; exit 1; }

mkdir -p $V2
REG=$V2/regulons.csv
AUC=$V2/auc_mtx.loom

# cisTarget DB: prefer v2, fall back to v1
pick() { for d in $V2/cistarget $V1/cistarget; do
  f=$(ls $d/$1 2>/dev/null | head -1); [[ -n "$f" ]] && { echo "$f"; return; }
done; }
FEATHER=$(pick '*.feather'); MOTIFS=$(pick '*motifs*.tbl')
[[ -n "$FEATHER" && -n "$MOTIFS" ]] || { echo "ERROR: cistarget DB missing"; exit 1; }
echo "Feather: $FEATHER"; echo "Motifs: $MOTIFS"; echo "ADJ: $ADJ"; echo "LOOM: $LOOM"

echo '=== ctx ==='
conda run -n epn2 pyscenic ctx "$ADJ" "$FEATHER" \
    --annotations_fname "$MOTIFS" --expression_mtx_fname "$LOOM" \
    --output "$REG" --mask_dropouts --num_workers 4

echo '=== aucell ==='
conda run -n epn2 pyscenic aucell "$LOOM" "$REG" \
    --output "$AUC" --num_workers 4

echo '=== downstream Fig 6 (v2) ==='
conda run -n epn2 python projectmd/v2/fig06_pyscenic_downstream.py

echo '✅ v2 pySCENIC pipeline + Fig 6 done'
