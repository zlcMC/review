#!/usr/bin/env python
"""fig06_regulon_subclone_differential.py — 新增 v2 分析（任务 B）.

针对 GTE009 Subclone_1 vs Subclone_2 做 regulon AUC 的 Mann-Whitney 差异检验，
把 Fig 1（subclone）与 Fig 6（regulon）打通。
也跨 4 个样本做 pooled Subclone_1 vs Subclone_2。

Inputs : output/fig6/auc_mtx.loom + scenic_input/meta.csv
Output : output/v2/fig6/
  - regulon_subclone_diff_GTE009.csv
  - regulon_subclone_diff_pooled.csv
  - Fig6_regulon_subclone_top_GTE009.pdf
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import dep_file, v2_dir, bh_adjust  # noqa

OUT = v2_dir("fig6")
auc_loom = Path(dep_file("fig6", "auc_mtx.loom"))
meta_csv = Path(dep_file("fig6", "scenic_input", "meta.csv"))


def decode(values):
    return [v.decode("utf-8") if isinstance(v, (bytes, bytearray)) else str(v) for v in values]


with h5py.File(str(auc_loom), "r") as h5:
    arr = h5["col_attrs"]["RegulonsAUC"][:]
    cells = decode(h5["col_attrs"]["CellID"][:])
auc = pd.DataFrame({n: arr[n] for n in arr.dtype.names}, index=cells, dtype=np.float32).T
print(f"AUC: {auc.shape[0]} regulons × {auc.shape[1]} cells")

meta = pd.read_csv(meta_csv, index_col=0)
common = auc.columns.intersection(meta.index)
auc, meta = auc.loc[:, common], meta.loc[common]


def regulon_diff(cells_g1, cells_g2, label1="Subclone_1", label2="Subclone_2"):
    a = auc.loc[:, cells_g1].to_numpy()
    b = auc.loc[:, cells_g2].to_numpy()
    n1, n2 = a.shape[1], b.shape[1]
    rows = []
    for i, reg in enumerate(auc.index):
        x, y = a[i], b[i]
        try:
            stat = mannwhitneyu(x, y, alternative="two-sided")
            p = float(stat.pvalue)
        except ValueError:
            p = np.nan
        rows.append({
            "regulon": reg,
            f"mean_{label1}": float(x.mean()),
            f"mean_{label2}": float(y.mean()),
            f"delta_{label1}_minus_{label2}": float(x.mean() - y.mean()),
            f"n_{label1}": n1,
            f"n_{label2}": n2,
            "mannwhitney_p": p,
        })
    df = pd.DataFrame(rows)
    df["fdr"] = bh_adjust(df["mannwhitney_p"].fillna(1.0).to_numpy())
    return df.sort_values("delta_" + label1 + "_minus_" + label2, ascending=False)


# --- GTE009 only ---
m9 = meta[(meta["sample"] == "GTE009") & meta["CNV_Cluster"].isin(["Subclone_1", "Subclone_2"])]
g1_cells = m9.index[m9["CNV_Cluster"] == "Subclone_1"]
g2_cells = m9.index[m9["CNV_Cluster"] == "Subclone_2"]
print(f"GTE009: Sub1={len(g1_cells)}, Sub2={len(g2_cells)}")
diff_g = regulon_diff(g1_cells, g2_cells)
diff_g.to_csv(OUT / "regulon_subclone_diff_GTE009.csv", index=False)

# --- pooled across 4 samples (only Subclone_1 vs Subclone_2) ---
mp = meta[meta["CNV_Cluster"].isin(["Subclone_1", "Subclone_2"])]
g1p = mp.index[mp["CNV_Cluster"] == "Subclone_1"]
g2p = mp.index[mp["CNV_Cluster"] == "Subclone_2"]
print(f"Pooled: Sub1={len(g1p)}, Sub2={len(g2p)}")
diff_p = regulon_diff(g1p, g2p)
diff_p.to_csv(OUT / "regulon_subclone_diff_pooled.csv", index=False)

# --- top regulons plot (GTE009) ---
top_up = diff_g.head(15)
top_dn = diff_g.tail(15).iloc[::-1]
plot_df = pd.concat([top_up, top_dn]).reset_index(drop=True)
plot_df["color"] = ["#b23a48" if d > 0 else "#3b6fa3"
                    for d in plot_df["delta_Subclone_1_minus_Subclone_2"]]
fig, ax = plt.subplots(figsize=(5.5, 8))
ax.barh(range(len(plot_df))[::-1], plot_df["delta_Subclone_1_minus_Subclone_2"],
        color=plot_df["color"], edgecolor="black", linewidth=0.4)
ax.set_yticks(range(len(plot_df))[::-1])
ax.set_yticklabels(plot_df["regulon"], fontsize=7)
ax.axvline(0, color="0.5", lw=0.6)
ax.set_xlabel("Mean AUC: Subclone_1 - Subclone_2")
ax.set_title("GTE009 differential regulons (top 15 up + 15 down)")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT / f"Fig6_regulon_subclone_top_GTE009.{ext}", dpi=300)
plt.close(fig)

n_sig_g = int((diff_g["fdr"] < 0.05).sum())
n_sig_p = int((diff_p["fdr"] < 0.05).sum())
print(f"GTE009 FDR<0.05: {n_sig_g}/{len(diff_g)}")
print(f"Pooled FDR<0.05: {n_sig_p}/{len(diff_p)}")
print(f"✓ → {OUT}/regulon_subclone_diff_*.csv & Fig6_regulon_subclone_top_GTE009.{{pdf,png}}")
