#!/bin/bash
# 30_infercnv_gte009.slurm — Fig 2D-F strict CNV candidate via inferCNV
#
# 用法:
#   cd /project/$USER/wslproject/projectmd/hpc
#   sbatch 30_infercnv_gte009.slurm

#SBATCH --job-name=infercnv_gte009
#SBATCH --partition=hpxg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=/project/%u/wslproject/logs/infercnv_gte009_%j.out
#SBATCH --error=/project/%u/wslproject/logs/infercnv_gte009_%j.err

set -euo pipefail

PROJECT_ROOT="/project/${USER}/wslproject"
GENE_ORDER="${PROJECT_ROOT}/output/reference/gencode_hg38_gene_order.tsv"
OUT_DIR="${PROJECT_ROOT}/output/fig2_infercnv_strict/GTE009"
RDS="${PROJECT_ROOT}/projectfile/GTE009_seurat.rds"

mkdir -p "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/output/reference" "${OUT_DIR}"

echo "节点: $(hostname)"
echo "开始时间: $(date)"
echo "核数: ${SLURM_CPUS_PER_TASK}"
echo "内存上限: ${SLURM_MEM_PER_NODE} MB"
free -h

source /project/${USER}/miniforge3/bin/activate
if conda env list | awk '{print $1}' | grep -qx infercnv_r; then
	conda activate infercnv_r
elif conda env list | awk '{print $1}' | grep -qx epn2_r; then
	conda activate epn2_r
else
	echo "ERROR: missing conda env infercnv_r or epn2_r. Run projectmd/hpc/01_prepare_infercnv_env.sh first."
	exit 1
fi

echo "R: $(which R)"
Rscript -e "cat(R.version.string, '\n'); cat('infercnv=', as.character(packageVersion('infercnv')), '\n'); cat('Seurat=', as.character(packageVersion('Seurat')), '\n')"

[[ -s "${RDS}" ]] || { echo "ERROR: 缺少 ${RDS}"; exit 1; }
[[ -s "${GENE_ORDER}" ]] || { echo "ERROR: 缺少 ${GENE_ORDER}"; echo "先用 fig02_infercnv_build_gene_order.py 从 GENCODE GTF 生成。"; exit 1; }

cd "${PROJECT_ROOT}"
Rscript projectmd/fig02_infercnv_run_gte009.R "${GENE_ORDER}" "${OUT_DIR}"

echo
echo "完成时间: $(date)"
find "${OUT_DIR}" -maxdepth 2 -type f -printf "%p\t%s bytes\n" | sort