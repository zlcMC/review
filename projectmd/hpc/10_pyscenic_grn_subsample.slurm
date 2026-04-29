#!/bin/bash
# 10_pyscenic_grn_subsample.slurm — pySCENIC GRN, 备用方案
# 当全量 35k 跑挂掉 (OOM 或 24h 超时) 时用此版本: 降采样到 15k 细胞
#
# 工作原理: 用 python 脚本从 integrated.loom 里抽 15k 细胞 → integrated_sub15k.loom,
#           然后跑 GRN, 输出到 fig6/adj_sub15k.tsv

#SBATCH --job-name=scenic_grn_sub
#SBATCH --partition=hpxg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=70G
#SBATCH --time=12:00:00
#SBATCH --output=/project/%u/wslproject/logs/scenic_grn_sub_%j.out
#SBATCH --error=/project/%u/wslproject/logs/scenic_grn_sub_%j.err

set -euo pipefail

PROJECT_ROOT="/project/${USER}/wslproject"
FIG6="${PROJECT_ROOT}/fig6"
LOOM_FULL="${FIG6}/scenic_input/integrated.loom"
LOOM_SUB="${FIG6}/scenic_input/integrated_sub15k.loom"
TFS="${FIG6}/cistarget/allTFs_hg38.txt"
ADJ="${FIG6}/adj_sub15k.tsv"
N_CELLS=15000

source /project/${USER}/miniforge3/bin/activate
conda activate scenic
echo "节点: $(hostname); 开始: $(date)"

# ---- 1) 降采样: 写一个新 loom ----
if [[ ! -s "${LOOM_SUB}" ]]; then
    echo "=== 降采样 → ${N_CELLS} 细胞 ==="
    python - <<EOF
import loompy, numpy as np
np.random.seed(42)
src = "${LOOM_FULL}"
dst = "${LOOM_SUB}"
n_target = ${N_CELLS}
with loompy.connect(src, mode='r') as ds:
    n_cells = ds.shape[1]
    idx = np.sort(np.random.choice(n_cells, size=min(n_target, n_cells), replace=False))
    print(f"原始: {ds.shape}; 采样到 {len(idx)} 细胞")
    mat = ds[:, idx]
    row_attrs = {k: ds.ra[k] for k in ds.ra.keys()}
    col_attrs = {k: ds.ca[k][idx] for k in ds.ca.keys()}
loompy.create(dst, mat, row_attrs, col_attrs)
print(f"✓ 写入 {dst}")
EOF
    ls -lh "${LOOM_SUB}"
fi

# ---- 2) 跑 GRN ----
if [[ -s "${ADJ}" ]]; then
    echo "✓ ADJ 已存在 (${ADJ}), 跳过 GRN"
    exit 0
fi

echo "=== GRN (subsample) ==="
pyscenic grn \
    "${LOOM_SUB}" "${TFS}" \
    --output "${ADJ}" \
    --num_workers ${SLURM_CPUS_PER_TASK} \
    --seed 42

echo
echo "✓ GRN (subsample) 完成: $(date)"
echo "ADJ 行数: $(wc -l < ${ADJ})"
ls -lh "${ADJ}"
echo
echo "提示: 用此 ADJ 跑 ctx/aucell 时, 把 20_pyscenic_ctx_aucell.slurm 里的"
echo "      ADJ 和 LOOM 路径改为 _sub15k 版本。"
