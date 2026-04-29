#!/usr/bin/env bash
# 00_prepare_conda_env.sh — 一次性: 在 /project 下装 miniforge + 建 scenic 环境
# 用法: 在 swarm01 登录节点上跑 `bash 00_prepare_conda_env.sh`
#
# 武大超算环境备注:
#   - 家目录 /home/$USER 只有 1 GB, 所有大文件必须放 /project
#   - 登录节点国际网络受限, 用清华镜像下 miniforge, 用清华源装 conda 包
#   - 装好后 `conda activate scenic` 永久可用
set -euo pipefail

PROJECT_ROOT="/project/${USER}"
MINIFORGE="${PROJECT_ROOT}/miniforge3"
ENV_NAME="scenic"

# ---- 1) 装 miniforge ----
if [[ ! -d "${MINIFORGE}" ]]; then
    cd "${PROJECT_ROOT}"
    echo "下载 miniforge (清华镜像)..."
    wget -q --show-progress \
        https://mirrors.tuna.tsinghua.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh \
        -O Miniforge3.sh
    bash Miniforge3.sh -b -p "${MINIFORGE}"
    rm Miniforge3.sh
    echo "✓ miniforge 装到 ${MINIFORGE}"
else
    echo "✓ miniforge 已存在: ${MINIFORGE}"
fi

# ---- 2) 配 conda 国内源 + 默认目录到 /project ----
mkdir -p "${PROJECT_ROOT}/conda_envs" "${PROJECT_ROOT}/conda_pkgs"
cat > ~/.condarc <<EOF
channels:
  - conda-forge
  - bioconda
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
show_channel_urls: true
envs_dirs:
  - ${PROJECT_ROOT}/conda_envs
pkgs_dirs:
  - ${PROJECT_ROOT}/conda_pkgs
EOF
echo "✓ ~/.condarc 已配置 (清华源 + envs_dirs 指向 /project)"

# ---- 3) 激活 miniforge ----
source "${MINIFORGE}/bin/activate"

# ---- 4) 建 scenic env ----
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "✓ conda env '${ENV_NAME}' 已存在, 跳过创建"
else
    echo "建 conda env '${ENV_NAME}' (pyscenic 0.12.1) ..."
    mamba create -n "${ENV_NAME}" -c conda-forge -c bioconda \
        python=3.10 pyscenic=0.12.1 \
        "dask<2023" "distributed<2023" \
        "setuptools<70" \
        numpy pandas scipy scikit-learn loompy pyarrow numba \
        -y
fi

# ---- 5) 测试 ----
conda activate "${ENV_NAME}"
echo
echo "=== 测试 ==="
echo "Python: $(which python)"
python --version
echo "pyscenic:"
pyscenic --help 2>&1 | head -10

echo
echo "✓ 环境就绪. 之后在 SLURM 脚本里用:"
echo "    source ${MINIFORGE}/bin/activate"
echo "    conda activate ${ENV_NAME}"
