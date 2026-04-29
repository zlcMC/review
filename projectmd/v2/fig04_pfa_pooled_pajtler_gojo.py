"""v2 add-on: pooled PF_EPN_A primary-vs-relapse trajectory meta
combining Pajtler GSE64415 bulk microarray + Gojo GSE141460 scRNA pseudobulk.

Outputs (output/v2/fig4_external_pajtler/):
  PFA_pooled_primary_vs_relapse_input.csv
  PFA_pooled_primary_vs_relapse_stats.csv
  Fig4_PFA_pooled_primary_vs_relapse.pdf/.png
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import v2_file, dep_file  # type: ignore

OUT_DIR = Path(v2_file("fig4_external_pajtler", "_dummy")).parent

# 1) Pajtler PF-A bulk: each sample is one row
pajtler = pd.read_csv(dep_file("fig4_external_pajtler", "GSE64415_trajectory_scores.csv"))
pajtler = pajtler[(pajtler["subgroup"] == "PF_EPN_A") &
                  (pajtler["condition"].isin(["primary", "relapse"]))].copy()
pajtler["cohort"] = "Pajtler_GSE64415_bulk"
pajtler["sample"] = pajtler["gsm"]
pajtler["condition_canon"] = pajtler["condition"].map({"primary": "Primary", "relapse": "Recurrent"})
pajtler_long = pajtler[["cohort", "sample", "condition_canon", "trajectory_score"]].rename(
    columns={"condition_canon": "condition"})

# 2) Gojo PF-A scRNA pseudobulk per sample (use median trajectory across malignant cells)
gojo_cells = pd.read_csv(dep_file(
    "fig4_gojo_scrna", "GSE141460_cell_metadata_with_trajectory_scores.csv"))
# subset to malignant + PF-A samples used in the prior PF-A vs recurrent fig
gojo_cells = gojo_cells[gojo_cells["malignant"].astype(str).str.lower().isin(["true", "1", "malignant", "yes"])].copy()
print(f"Gojo malignant cells: {len(gojo_cells)}")

gojo_pfa = gojo_cells[gojo_cells["subtype"].astype(str).str.upper().isin(["PFA", "PF-A", "PF_A"])].copy()
print(f"Gojo PF-A malignant cells: {len(gojo_pfa)}; samples: {gojo_pfa['sample'].nunique()}")
gojo_per_sample = gojo_pfa.groupby(["sample", "condition"], dropna=True)["trajectory_score"].median().reset_index()
gojo_per_sample["cohort"] = "Gojo_GSE141460_scRNA"
gojo_per_sample["condition"] = gojo_per_sample["condition"].astype(str).str.title()
gojo_long = gojo_per_sample[["cohort", "sample", "condition", "trajectory_score"]]

# 3) Pool & write
pool = pd.concat([pajtler_long, gojo_long], ignore_index=True)
pool = pool[pool["condition"].isin(["Primary", "Recurrent"])].copy()
pool.to_csv(OUT_DIR / "PFA_pooled_primary_vs_relapse_input.csv", index=False)
print("\nPooled PF-A samples per cohort/condition:")
print(pool.groupby(["cohort", "condition"]).size())

# 4) Stats
def wilcox(a, b):
    if len(a) > 0 and len(b) > 0:
        try:
            return stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError:
            return float("nan")
    return float("nan")

p_overall = wilcox(pool.loc[pool.condition=="Primary", "trajectory_score"],
                   pool.loc[pool.condition=="Recurrent", "trajectory_score"])

# cohort-adjusted linear model: trajectory ~ recurrent + cohort
df = pool.copy()
df["recurrent"] = (df.condition == "Recurrent").astype(int)
df_dum = pd.get_dummies(df[["recurrent", "cohort"]], columns=["cohort"], drop_first=True, dtype=float)
X = sm.add_constant(df_dum)
y = df["trajectory_score"].astype(float)
model = sm.OLS(y, X).fit()
print(model.summary())

stats_rows = [
    {"test": "wilcoxon_pooled_unadj",
     "p_value": p_overall,
     "n_primary": int((pool.condition=="Primary").sum()),
     "n_recurrent": int((pool.condition=="Recurrent").sum())},
    {"test": "lm_cohort_adjusted",
     "p_value": float(model.pvalues["recurrent"]),
     "estimate": float(model.params["recurrent"]),
     "n_primary": int((pool.condition=="Primary").sum()),
     "n_recurrent": int((pool.condition=="Recurrent").sum())},
]
pd.DataFrame(stats_rows).to_csv(OUT_DIR / "PFA_pooled_primary_vs_relapse_stats.csv", index=False)

# 5) Plot
fig, ax = plt.subplots(figsize=(6.0, 4.6))
sns.boxplot(data=pool, x="condition", y="trajectory_score",
            order=["Primary", "Recurrent"], color="white", width=0.45, fliersize=0, ax=ax)
sns.stripplot(data=pool, x="condition", y="trajectory_score",
              order=["Primary", "Recurrent"], hue="cohort",
              palette={"Pajtler_GSE64415_bulk": "#2C7FB8", "Gojo_GSE141460_scRNA": "#E08214"},
              size=5, jitter=0.18, ax=ax)
ax.set_title(
    f"PF_EPN_A primary vs recurrent (Pajtler bulk + Gojo scRNA)\n"
    f"Wilcoxon p={p_overall:.3g} | cohort-adjusted lm p={model.pvalues['recurrent']:.3g}",
    fontsize=11)
ax.set_xlabel("")
ax.set_ylabel("Trajectory score (undiff - diff)")
ax.legend(loc="best", fontsize=8, frameon=True)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT_DIR / f"Fig4_PFA_pooled_primary_vs_relapse.{ext}", dpi=300)
plt.close(fig)
print("\nDone. Outputs in:", OUT_DIR)
