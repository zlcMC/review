#!/usr/bin/env python3
"""task 3 — subtype-matched 跨队列 DEG/LR concordance（拆开 confounder 报告）.

策略：
- Gillen GSE125969 recurrent 全是 ST-RELA → 用 ST-RELA-only DEG
- Gojo GSE141460 PFA-only DEG 已存在
- 双向 concordance:
  (a) Gillen ST-RELA-only vs Gojo ST-RELA-only（如 Gojo 有 ST-RELA primary 与 recurrent）
  (b) Gillen ST-RELA-only vs Gojo PFA-only（cross-subtype，预期低）
  (c) Gojo PFA-only vs Wu_GTE009 Subclone1vs2 DEG（不同维度，作为 sanity）
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import dep_file, v2_dir  # noqa

OUT = v2_dir("fig4_external_scrna_meta")


def load_deg(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    delta = "recurrent_minus_primary_logFC" if "recurrent_minus_primary_logFC" in df.columns else "recurrent_minus_primary"
    return pd.DataFrame({"gene": df["gene"].astype(str), f"{label}_delta": df[delta]})


def concordance(a: pd.DataFrame, b: pd.DataFrame, label_a: str, label_b: str) -> dict:
    m = a.merge(b, on="gene", how="inner")
    m["sa"] = np.sign(m[f"{label_a}_delta"])
    m["sb"] = np.sign(m[f"{label_b}_delta"])
    nz = m[(m["sa"] != 0) & (m["sb"] != 0)].copy()
    nz["concordant"] = nz["sa"] == nz["sb"]
    n_conc = int(nz["concordant"].sum())
    n_tot = len(nz)
    p = float(binomtest(n_conc, n_tot, 0.5, alternative="greater").pvalue) if n_tot else np.nan
    rho, rho_p = (np.nan, np.nan) if n_tot < 5 else spearmanr(nz[f"{label_a}_delta"], nz[f"{label_b}_delta"])
    return {
        "comparison": f"{label_a} vs {label_b}",
        "n_genes_intersect": n_tot,
        "concordant": n_conc,
        "fraction_concordant": n_conc / n_tot if n_tot else np.nan,
        "binomial_p_greater_than_half": p,
        "spearman_rho": float(rho) if not np.isnan(rho) else np.nan,
        "spearman_p": float(rho_p) if not np.isnan(rho_p) else np.nan,
    }


# Load DEG sets
gillen_strela = load_deg(Path(dep_file("fig4_scrna", "GSE125969_neoplastic_pseudobulk_DEG_ST_RELA_only.csv")), "Gillen_STRELA")
gillen_all    = load_deg(Path(dep_file("fig4_scrna", "GSE125969_neoplastic_pseudobulk_DEG_all_subtypes.csv")), "Gillen_all")
gojo_pfa      = load_deg(Path(dep_file("fig4_gojo_scrna", "GSE141460_malignant_pseudobulk_DEG_PFA_only.csv")), "Gojo_PFA")
gojo_all      = load_deg(Path(dep_file("fig4_gojo_scrna", "GSE141460_malignant_pseudobulk_DEG_all_subtypes.csv")), "Gojo_all")

results = [
    # 同 subtype 跨队列：受限于双方各自只有一个 subtype
    concordance(gillen_strela, gojo_all,    "Gillen_STRELA", "Gojo_all"),
    concordance(gillen_all,    gojo_pfa,    "Gillen_all",    "Gojo_PFA"),
    # 跨 subtype（参照）：恰好两方各自只剩自己 subtype
    concordance(gillen_strela, gojo_pfa,    "Gillen_STRELA", "Gojo_PFA"),
    # baseline pooled（已有，重算放一起）
    concordance(gillen_all,    gojo_all,    "Gillen_all",    "Gojo_all"),
]
df = pd.DataFrame(results)
df.to_csv(OUT / "subtype_matched_DEG_concordance.csv", index=False)
print(df.to_string(index=False))
print(f"\n→ {OUT}/subtype_matched_DEG_concordance.csv")
