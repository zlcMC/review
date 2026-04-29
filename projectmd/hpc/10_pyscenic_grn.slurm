#!/bin/bash
# 10_pyscenic_grn.slurm — pySCENIC step 1: GRN (GRNBoost2)
# 全量 35k 细胞 × 16k 基因, 内存压力大
#
# 用法: cd /project/$USER/wslproject/projectmd/hpc && sbatch 10_pyscenic_grn.slurm

#SBATCH --job-name=scenic_grn
#SBATCH --partition=hpxg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=/project/%u/wslproject/logs/scenic_grn_%j.out
#SBATCH --error=/project/%u/wslproject/logs/scenic_grn_%j.err

set -euo pipefail

PROJECT_ROOT="/project/${USER}/wslproject"
FIG6="${PROJECT_ROOT}/fig6"
LOOM="${FIG6}/scenic_input/integrated.loom"
TFS="${FIG6}/cistarget/allTFs_hg38.txt"
ADJ="${FIG6}/adj.tsv"

# ---- 环境: 激活 miniforge 里的 scenic env ----
source /project/${USER}/miniforge3/bin/activate
conda activate scenic
echo "节点: $(hostname)"
echo "Python: $(which python)"
echo "pyscenic: $(pyscenic --version 2>&1 | head -1)"
echo "核数: ${SLURM_CPUS_PER_TASK}"
echo "内存上限: ${SLURM_MEM_PER_NODE} MB"
echo "开始时间: $(date)"
free -h

# ---- 检查输入 ----
[[ -s "${LOOM}" ]] || { echo "ERROR: 缺少 ${LOOM}"; exit 1; }
[[ -s "${TFS}"  ]] || { echo "ERROR: 缺少 ${TFS}"; exit 1; }

if [[ -s "${ADJ}" ]]; then
    echo "✓ ADJ 已存在 (${ADJ}), 跳过 GRN"
    exit 0
fi

# ---- 跑 GRN ----
pyscenic grn \
    "${LOOM}" "${TFS}" \
    --output "${ADJ}" \
    --num_workers ${SLURM_CPUS_PER_TASK} \
    --seed 42

echo
echo "✓ GRN 完成: $(date)"
echo "ADJ 行数: $(wc -l < ${ADJ})"
ls -lh "${ADJ}"
