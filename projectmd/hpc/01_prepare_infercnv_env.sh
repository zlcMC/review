#!/bin/bash
# Prepare a lightweight R environment for the GTE009 inferCNV job on WHU HPC.

set -euo pipefail

ENV_NAME="infercnv_r"
MINIFORGE="/project/${USER}/miniforge3"

if [[ ! -s "${MINIFORGE}/bin/activate" ]]; then
    echo "ERROR: Miniforge not found at ${MINIFORGE}. Run 00_prepare_conda_env.sh first."
    exit 1
fi

source "${MINIFORGE}/bin/activate"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "✓ conda env ${ENV_NAME} already exists"
else
    if command -v mamba >/dev/null 2>&1; then
        SOLVER=mamba
    else
        SOLVER=conda
    fi

    echo "Creating ${ENV_NAME} with ${SOLVER}..."
    "${SOLVER}" create -y -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        'r-base>=4.3,<4.5' \
        r-seurat \
        r-seuratobject \
        r-hdf5r \
        r-data.table \
        bioconductor-infercnv
fi

conda activate "${ENV_NAME}"
Rscript -e "cat(R.version.string, '\n'); cat('Seurat=', as.character(packageVersion('Seurat')), '\n'); cat('infercnv=', as.character(packageVersion('infercnv')), '\n')"

echo "✓ ${ENV_NAME} is ready"