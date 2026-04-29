"""v2 final integration: 4-cohort PF-A primary vs relapse trajectory panel.

Combines:
- Wu (this paper) GSE189289: 4 PF-A primary samples (cell-level scores in
  output/fig3/GTE*_scores.rds; aggregate to per-sample median)
- Gillen GSE125969: PF-A subset from
  output/v2/fig4_scrna/GSE125969_neoplastic_cell_trajectory_scores.csv
- Gojo GSE141460: PF-A subset from
  output/v2/fig4_gojo_scrna/GSE141460_malignant_cell_trajectory_scores.csv
- Pajtler GSE64415 bulk: PF-A from
  output/v2/fig4_external_pajtler/GSE64415_trajectory_scores.csv

Outputs (output/v2/fig4_external_pajtler/):
  Fig4D_final_4cohort_PFA_trajectory_per_sample.csv
  Fig4D_final_4cohort_PFA_trajectory.pdf/.png
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
from helpers.helpers_common import v2_file, dep_file, output_path  # type: ignore

OUT_DIR = Path(v2_file("fig4_external_pajtler", "_dummy")).parent

# ---- 1) Wu: read 4 R rds via subprocess (call Rscript to dump CSV)
WU_TMP = OUT_DIR / "wu_per_cell_trajectory.csv"
import subprocess
r_code = f"""
samples <- c('GTE001','GTE002','GTE009','GTE012')
out <- do.call(rbind, lapply(samples, function(s) {{
  rds <- file.path('output','fig3', paste0(s, '_scores.rds'))
  if (!file.exists(rds)) return(NULL)
  d <- readRDS(rds)
  data.frame(sample=s, trajectory_score=d$trajectory_score)
}}))
write.csv(out, '{WU_TMP}', row.names=FALSE)
cat('Wu rows:', nrow(out), '\\n')
"""
subprocess.run(["conda", "run", "-n", "epn2_r", "Rscript", "-e", r_code], check=True)
wu_cells = pd.read_csv(WU_TMP)
wu_per_sample = wu_cells.groupby("sample")["trajectory_score"].median().reset_index()
wu_per_sample["cohort"] = "Wu_GSE189289_scRNA"
wu_per_sample["condition"] = "Primary"  # all 4 are PF-A primary
wu_per_sample["subtype"] = "PF-A"
print(f"Wu per-sample (PF-A primary): {len(wu_per_sample)}")

# ---- 2) Gillen
g_neo = pd.read_csv(dep_file("fig4_scrna", "GSE125969_neoplastic_cell_trajectory_scores.csv"))
print("Gillen neo cols:", list(g_neo.columns)[:15])
# determine subtype col
sub_col_g = next((c for c in g_neo.columns if "subtype" in c.lower() or "subgroup" in c.lower() or "molecular" in c.lower()), None)
cond_col_g = next((c for c in g_neo.columns if c.lower() in ("condition", "primary_recurrent", "status", "presentation")), None)
samp_col_g = next((c for c in g_neo.columns if c.lower() in ("sample", "sample_id", "patient", "patient_id")), None)
print(f"Gillen picks: subtype={sub_col_g} cond={cond_col_g} sample={samp_col_g}")
if sub_col_g is None:
    g_pfa = pd.DataFrame()
else:
    g_pfa = g_neo[g_neo[sub_col_g].astype(str).str.upper().str.contains("PFA|PF.A|PF_A", regex=True, na=False)].copy()
    # Gillen presentation: 1=Primary, 2=Recurrent (per Gillen 2020 readme)
    g_pfa["condition"] = g_pfa[cond_col_g].map({1: "Primary", 2: "Recurrent"})
    g_per_sample = g_pfa.groupby([samp_col_g, "condition"], dropna=True)["trajectory_score"].median().reset_index()
    g_per_sample = g_per_sample.rename(columns={samp_col_g: "sample"})
    g_per_sample["cohort"] = "Gillen_GSE125969_scRNA"
    g_per_sample["subtype"] = "PF-A"
print(f"Gillen PF-A per-sample: {len(g_per_sample) if sub_col_g else 0}")

# ---- 3) Gojo
gojo_cells = pd.read_csv(dep_file("fig4_gojo_scrna", "GSE141460_cell_metadata_with_trajectory_scores.csv"))
gojo_mal = gojo_cells[gojo_cells["malignant"].astype(str).str.lower().isin(["true","1","malignant","yes"])]
gojo_pfa = gojo_mal[gojo_mal["subtype"].astype(str).str.upper().isin(["PFA","PF-A","PF_A"])].copy()
gojo_per_sample = gojo_pfa.groupby(["sample", "condition"])["trajectory_score"].median().reset_index()
gojo_per_sample["cohort"] = "Gojo_GSE141460_scRNA"
gojo_per_sample["subtype"] = "PF-A"
gojo_per_sample["condition"] = gojo_per_sample["condition"].astype(str).str.title()
print(f"Gojo PF-A per-sample: {len(gojo_per_sample)}")

# ---- 4) Pajtler
paj = pd.read_csv(dep_file("fig4_external_pajtler", "GSE64415_trajectory_scores.csv"))
paj_pfa = paj[(paj["subgroup"]=="PF_EPN_A") & (paj["condition"].isin(["primary","relapse"]))].copy()
paj_pfa = paj_pfa.rename(columns={"gsm": "sample"})
paj_pfa["cohort"] = "Pajtler_GSE64415_bulk"
paj_pfa["subtype"] = "PF-A"
paj_pfa["condition"] = paj_pfa["condition"].map({"primary":"Primary", "relapse":"Recurrent"})
paj_pfa = paj_pfa[["cohort", "sample", "condition", "subtype", "trajectory_score"]]
print(f"Pajtler PF-A per-sample: {len(paj_pfa)}")

# ---- 5) Pool
pool = pd.concat([
    wu_per_sample[["cohort","sample","condition","subtype","trajectory_score"]],
    g_per_sample[["cohort","sample","condition","subtype","trajectory_score"]] if sub_col_g else pd.DataFrame(),
    gojo_per_sample[["cohort","sample","condition","subtype","trajectory_score"]],
    paj_pfa,
], ignore_index=True)
pool = pool[pool["condition"].isin(["Primary", "Recurrent"])].copy()
pool.to_csv(OUT_DIR / "Fig4D_final_4cohort_PFA_trajectory_per_sample.csv", index=False)

print("\nPer cohort × condition counts:")
print(pool.groupby(["cohort", "condition"]).size())

# ---- 6) Stats: pooled cohort-adjusted lm
df = pool.copy()
df["recurrent"] = (df.condition=="Recurrent").astype(int)
df_dum = pd.get_dummies(df[["recurrent","cohort"]], columns=["cohort"], drop_first=True, dtype=float)
X = sm.add_constant(df_dum)
y = df["trajectory_score"].astype(float)
m = sm.OLS(y, X).fit()
print(m.summary())
n_p = int((pool.condition=="Primary").sum())
n_r = int((pool.condition=="Recurrent").sum())

# ---- 7) Plot
plt.rcParams["font.family"] = "DejaVu Sans"
fig, ax = plt.subplots(figsize=(7.6, 5.0))
order = ["Primary", "Recurrent"]
cohort_order = ["Wu_GSE189289_scRNA", "Gillen_GSE125969_scRNA", "Gojo_GSE141460_scRNA", "Pajtler_GSE64415_bulk"]
palette = {
    "Wu_GSE189289_scRNA": "#222222",
    "Gillen_GSE125969_scRNA": "#66A61E",
    "Gojo_GSE141460_scRNA": "#E08214",
    "Pajtler_GSE64415_bulk": "#7570B3",
}
sns.boxplot(data=pool, x="condition", y="trajectory_score", order=order,
            color="white", width=0.5, fliersize=0, ax=ax)
sns.stripplot(data=pool, x="condition", y="trajectory_score", order=order,
              hue="cohort", hue_order=cohort_order, palette=palette,
              size=6, jitter=0.2, ax=ax, alpha=0.85, edgecolor="black", linewidth=0.5)
ax.set_xlabel("")
ax.set_ylabel("Trajectory score per sample (median)")
ax.set_title(
    f"Fig 4D reproduction: PF-A primary vs recurrent across 4 cohorts\n"
    f"n={n_p} primary vs {n_r} recurrent | cohort-adjusted lm p={m.pvalues['recurrent']:.3g}, β={m.params['recurrent']:+.3f}",
    fontsize=11)
ax.legend(title="Cohort", loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(OUT_DIR / f"Fig4D_final_4cohort_PFA_trajectory.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"\nDone. Saved to {OUT_DIR}")
