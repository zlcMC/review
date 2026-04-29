#!/usr/bin/env python
"""v2 of fig06_pyscenic_downstream.py — Fig 6 regulon AUC heatmap + cell-type specific regulons.

Input : output/fig6/auc_mtx.loom (from fig06_pyscenic_run.py — HPC, frozen)
Output: output/v2/fig6/Fig6_regulon_AUC_heatmap.pdf, regulon_celltype_top10.csv, ...
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import dep_file, v2_dir

OUT = v2_dir("fig6")
auc_loom = Path(dep_file("fig6", "auc_mtx.loom"))
meta_csv = Path(dep_file("fig6", "scenic_input", "meta.csv"))

if not auc_loom.exists():
    sys.exit(f"missing {auc_loom}; run fig06_pyscenic_run.py on the HPC first")


def decode(values):
    return [v.decode("utf-8") if isinstance(v, (bytes, bytearray)) else str(v) for v in values]


# pySCENIC's aucell loom stores regulon activities under col_attrs/RegulonsAUC;
# the main matrix is still the gene × cell expression matrix.
with h5py.File(str(auc_loom), "r") as h5:
    if "RegulonsAUC" not in h5["col_attrs"]:
        sys.exit(f"missing col_attrs/RegulonsAUC in {auc_loom}")
    regulons_auc = h5["col_attrs"]["RegulonsAUC"][:]
    regulon_names = list(regulons_auc.dtype.names or [])
    cell_ids = decode(h5["col_attrs"]["CellID"][:])
    if not regulon_names:
        sys.exit(f"empty RegulonsAUC in {auc_loom}")
    auc = pd.DataFrame(
        {name: regulons_auc[name] for name in regulon_names},
        index=cell_ids,
        dtype=np.float32,
    ).T
print(f"AUC: {auc.shape}")

meta = pd.read_csv(meta_csv, index_col=0)
common = auc.columns.intersection(meta.index)
auc = auc.loc[:, common]
meta = meta.loc[common]
ct_col = next(
    (c for c in ["Brief_cluster", "celltype", "cluster"] if c in meta.columns),
    meta.columns[0],
)
print(f"using {ct_col} ({meta[ct_col].nunique()} groups)")

mean_auc = auc.T.groupby(meta[ct_col]).mean().T  # regulon × celltype
mean_auc.to_csv(OUT / "regulon_AUC_mean_by_celltype.csv")

z = mean_auc.sub(mean_auc.mean(axis=1), axis=0).div(mean_auc.std(axis=1) + 1e-6, axis=0)
top = {ct: z.nlargest(10, ct).index.tolist() for ct in z.columns}
pd.DataFrame.from_dict(top, orient="index").T.to_csv(OUT / "regulon_celltype_top10.csv")

top_regs = sorted({g for v in top.values() for g in v})
plt.figure(figsize=(max(6, mean_auc.shape[1] * 0.5), max(8, len(top_regs) * 0.18)))
sns.heatmap(z.loc[top_regs], cmap="RdBu_r", center=0, cbar_kws={"label": "Z-score AUC"})
plt.title("Cell-type specific regulons (SCENIC AUC)")
plt.tight_layout()
plt.savefig(OUT / "Fig6_regulon_AUC_heatmap.pdf")
plt.close()
print(f"✓ Fig 6 → {OUT}/Fig6_regulon_AUC_heatmap.pdf")
