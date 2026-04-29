#!/usr/bin/env python
"""task 5 — cohort × subtype × condition 透明化交叉表（Fig4 supp table）.

汇总所有外部 + 本队列样本，以 cohort/subtype/condition 三维交叉表呈现，
帮助解释 Fig4 外部验证为何受 confounder 限制。
"""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import dep_file, v2_dir, subtype_group  # noqa

OUT = v2_dir("fig4_external_scrna_meta")

rows = []

# Wu 2022 GSE189289 本队列
for s in ["GTE001", "GTE002", "GTE009", "GTE012"]:
    rows.append({"cohort": "Wu_GSE189289", "sample_id": s,
                 "subtype": "PF-A", "subtype_group": "PF-A", "condition": "Primary"})

# Gillen GSE125969 scRNA
g = pd.read_csv(dep_file("fig4_scrna", "GSE125969_patient_trajectory_summary.csv"))
for _, r in g.iterrows():
    rows.append({"cohort": "Gillen_GSE125969", "sample_id": str(r["patient_id"]),
                 "subtype": r["tumor_subtype"], "subtype_group": subtype_group(r["tumor_subtype"]),
                 "condition": r["condition"]})

# Gojo GSE141460 scRNA
gj = pd.read_csv(dep_file("fig4_gojo_scrna", "GSE141460_sample_trajectory_summary.csv"))
for _, r in gj.iterrows():
    rows.append({"cohort": "Gojo_GSE141460", "sample_id": r["sample_file"],
                 "subtype": r["subtype"], "subtype_group": subtype_group(r["subtype"]),
                 "condition": r["condition"]})

df = pd.DataFrame(rows)
df.to_csv(OUT / "supp_cohort_subtype_condition_samples.csv", index=False)

# 三维交叉表
ct = pd.crosstab([df["cohort"], df["condition"]], df["subtype_group"], margins=True)
ct.to_csv(OUT / "supp_cohort_subtype_condition_crosstab.csv")
print("=== cohort × condition × subtype_group ===")
print(ct.to_string())
print(f"\n→ {OUT}/supp_cohort_subtype_condition_*.csv")
