#!/usr/bin/env python3
"""task 6 — chr11 4-sample expression-level validation.

不重跑 inferCNV (HPC only)，而是用 GTE001/002/009/012 的 h5ad 表达层
比较 chr11 基因（DYNC2H1 + 其它 chr11 cilium 基因）在 Subclone_1 vs Subclone_2 的
mean log-normalized expression delta。
回答：DYNC2H1 / chr11 在 Subclone_1 高的方向是否在 4 个样本里一致？
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import raw_data_path, v2_dir, dep_file  # noqa

OUT = v2_dir("fig2_infercnv_strict")
OUT.mkdir(parents=True, exist_ok=True)

# 用 v1 inferCNV 的 chr11 stats 拿到 chr11 基因清单（已存在）
chr11_csv = Path(dep_file("fig2_infercnv_strict", "GTE009_downstream", "GTE009_chr11_inferCNV_genelevel_stats.csv"))
chr11_genes = pd.read_csv(chr11_csv)["gene"].astype(str).tolist()
print(f"chr11 genes from inferCNV table: {len(chr11_genes)}")

samples = ["GTE001", "GTE002", "GTE009", "GTE012"]
records = []
for s in samples:
    h5 = raw_data_path(f"{s}.h5ad")
    if not Path(h5).exists():
        print(f"[skip] {h5} not found"); continue
    a = sc.read_h5ad(h5)
    # 提取 subclone 标签（CNV_Cluster 或 subclone 列）
    # 优先 'CNV_Cluster'（值是 Subclone_*），不是 'CNV_cluster'（值是数字 / Undefined）
    if "CNV_Cluster" in a.obs.columns:
        col = "CNV_Cluster"
    else:
        cand = [c for c in a.obs.columns if "subclone" in c.lower()]
        if not cand:
            print(f"[skip] {s}: no subclone column"); continue
        col = cand[0]
    a.obs[col] = a.obs[col].astype(str)
    # 仅 Subclone_1 / Subclone_2
    keep = a.obs[col].isin(["Subclone_1", "Subclone_2"])
    if keep.sum() < 50:
        print(f"[skip] {s}: too few subclone cells"); continue
    aa = a[keep].copy()
    # use raw counts → log-normalize (per current adata.X assumed log1p; if not, normalize)
    if aa.X.max() > 50:  # likely raw counts
        sc.pp.normalize_total(aa, target_sum=1e4)
        sc.pp.log1p(aa)
    common = [g for g in chr11_genes if g in aa.var_names]
    g1 = aa[aa.obs[col] == "Subclone_1"]
    g2 = aa[aa.obs[col] == "Subclone_2"]
    if g1.n_obs < 10 or g2.n_obs < 10:
        print(f"[skip] {s}: subclone size too small"); continue
    X1 = np.asarray(g1[:, common].X.mean(axis=0)).ravel()
    X2 = np.asarray(g2[:, common].X.mean(axis=0)).ravel()
    df = pd.DataFrame({
        "sample": s, "gene": common,
        "mean_Subclone_1": X1, "mean_Subclone_2": X2,
        "delta": X1 - X2,
        "n_Sub1": g1.n_obs, "n_Sub2": g2.n_obs,
    })
    records.append(df)
    print(f"  {s}: {len(common)}/{len(chr11_genes)} chr11 genes; Sub1={g1.n_obs} Sub2={g2.n_obs}")

if not records:
    sys.exit("No samples processed")

allrec = pd.concat(records, ignore_index=True)
allrec.to_csv(OUT / "chr11_expression_subclone_4samples_long.csv", index=False)

wide = allrec.pivot_table(index="gene", columns="sample", values="delta")
wide["sign_concordant_count"] = (wide.fillna(0) > 0).sum(axis=1)
wide["sign_dn_count"] = (wide.fillna(0) < 0).sum(axis=1)
wide.to_csv(OUT / "chr11_expression_subclone_4samples_wide.csv")

# DYNC2H1 single-row summary
print("\n=== DYNC2H1 across samples (delta = Sub1 - Sub2 mean log-norm expr) ===")
print(allrec[allrec["gene"] == "DYNC2H1"][["sample", "mean_Subclone_1", "mean_Subclone_2", "delta"]].to_string(index=False))

# Direction concordance per gene across samples (vote)
n_samp = wide.notna().sum(axis=1)
up_frac = wide.iloc[:, :len(samples)].apply(lambda r: (r > 0).sum() / r.notna().sum() if r.notna().sum() > 0 else np.nan, axis=1)
summary = pd.DataFrame({
    "gene": wide.index,
    "n_samples_with_data": n_samp.values,
    "frac_samples_Sub1_higher": up_frac.values,
})
summary.to_csv(OUT / "chr11_expression_subclone_4samples_summary.csv", index=False)

n_total = (summary["n_samples_with_data"] >= 3).sum()
n_sub1_consensus = ((summary["n_samples_with_data"] >= 3) & (summary["frac_samples_Sub1_higher"] >= 0.75)).sum()
n_sub2_consensus = ((summary["n_samples_with_data"] >= 3) & (summary["frac_samples_Sub1_higher"] <= 0.25)).sum()
print(f"\n=== chr11 gene-level consensus ===")
print(f"Total genes with data in >=3 samples: {n_total}")
print(f"  Sub1-higher in >=75% samples: {n_sub1_consensus}")
print(f"  Sub2-higher in >=75% samples: {n_sub2_consensus}")
print(f"\n→ {OUT}/chr11_expression_subclone_4samples_*.csv")
