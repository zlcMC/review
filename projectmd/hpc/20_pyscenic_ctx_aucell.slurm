#!/bin/bash
# 20_pyscenic_ctx_aucell.slurm — pySCENIC step 2+3: ctx + aucell
# 注意: ctx 的内存峰值取决于 num_workers (每个 worker 加载完整 ranking db ~3GB)
# 历史教训: 16 worker + 60G 在 1.5M 行 adj.tsv 上 OOM + 6h timeout
# 现配置: 8 worker + 80G + 24h, 速度损失 ~30% 但稳定完成
#
# 依赖: 需要 fig6/adj.tsv 已生成 (先跑 10_pyscenic_grn.slurm)
# 用法: sbatch 20_pyscenic_ctx_aucell.slurm

#SBATCH --job-name=scenic_ctx
#SBATCH --partition=hpxg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=/project/%u/wslproject/logs/scenic_ctx_%j.out
#SBATCH --error=/project/%u/wslproject/logs/scenic_ctx_%j.err

set -euo pipefail

PROJECT_ROOT="/project/${USER}/wslproject"
FIG6="${PROJECT_ROOT}/fig6"
LOOM="${FIG6}/scenic_input/integrated.loom"
ADJ="${FIG6}/adj.tsv"
RANK="${FIG6}/cistarget/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather"
MOTIFS="${FIG6}/cistarget/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
REG="${FIG6}/regulons.csv"
AUC="${FIG6}/auc_mtx.loom"

# ---- 环境 ----
source /project/${USER}/miniforge3/bin/activate
conda activate scenic
echo "节点: $(hostname); 开始: $(date)"
echo "pyscenic: $(pyscenic --version 2>&1 | head -1)"

[[ -s "${ADJ}"    ]] || { echo "ERROR: 缺少 ${ADJ}, 先跑 GRN"; exit 1; }
[[ -s "${RANK}"   ]] || { echo "ERROR: 缺少 ${RANK}";          exit 1; }
[[ -s "${MOTIFS}" ]] || { echo "ERROR: 缺少 ${MOTIFS}";        exit 1; }
[[ -s "${LOOM}"   ]] || { echo "ERROR: 缺少 ${LOOM}";          exit 1; }

# ctx 用少 worker 节省内存; aucell 排序统计无此问题, 用满核
CTX_WORKERS=8
AUCELL_WORKERS=${SLURM_CPUS_PER_TASK}

# ---- ctx ----
# 删掉空 regulons.csv (上次失败留下的 0 字节文件), 否则会被误认为已完成
[[ -f "${REG}" && ! -s "${REG}" ]] && rm -f "${REG}"

if [[ ! -s "${REG}" ]]; then
    echo "=== ctx (workers=${CTX_WORKERS}) ==="
    pyscenic ctx \
        "${ADJ}" "${RANK}" \
        --annotations_fname "${MOTIFS}" \
        --expression_mtx_fname "${LOOM}" \
        --output "${REG}" \
        --mask_dropouts \
        --num_workers ${CTX_WORKERS}
    echo "✓ regulons: $(wc -l < ${REG}) lines"
else
    echo "✓ regulons.csv 已存在, 跳过 ctx"
fi

# ---- aucell ----
if [[ ! -s "${AUC}" ]]; then
    echo "=== aucell (workers=${AUCELL_WORKERS}) ==="
    pyscenic aucell \
        "${LOOM}" "${REG}" \
        --output "${AUC}" \
        --num_workers ${AUCELL_WORKERS}
else
    echo "✓ auc_mtx.loom 已存在, 跳过 aucell"
fi

echo
echo "✓ 全部完成: $(date)"
ls -lh "${REG}" "${AUC}"
echo
echo "下载结果到本地 WSL:"
echo "  rsync -avzP zhanglichuan@swarm-xfe.whu.edu.cn:${REG} ${AUC} /home/zlcmc/wslproject/output/fig6/"
