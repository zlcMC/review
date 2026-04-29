#!/usr/bin/env python3
"""Gojo GSE141460 external scRNA validation for trajectory and LR panels."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, ttest_ind

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import (  # noqa: E402
    KEY_LR_GENES_UNION,
    OUTPUT_DIR,
    bh_adjust,
    dep_file,
    parse_tsv_floats,
    v2_dir,
    v2_file,
)


CACHE_DIR = OUTPUT_DIR / "external_cache" / "GSE141460"   # large cache stays at v1 location
ROOT_DIR = CACHE_DIR / "tpm_meta" / "Ependymoma"
PMC_DIR = OUTPUT_DIR / "external_cache" / "pmc_supp"     # large cache stays at v1 location
OUT_DIR = v2_dir("fig4_gojo_scrna")

KEY_LR_GENES = list(KEY_LR_GENES_UNION)


def condition_from_text(value: object) -> str:
    text = str(value).lower()
    if "diagnostic" in text:
        return "Primary"
    if "recurrence" in text or "progression" in text or "metastasis" in text:
        return "Recurrent"
    return "Unknown"


def condition_from_sample(sample: str) -> str | None:
    lower = sample.lower()
    if "dia" in lower:
        return "Primary"
    if "rec" in lower:
        return "Recurrent"
    if lower.startswith("muv043") or lower == "peds4":
        return "Recurrent"
    return None


def load_gojo_clinical() -> pd.DataFrame:
    path = PMC_DIR / "Gojo_mmc2.xlsx"
    if not path.exists():
        raise FileNotFoundError(f"Missing Gojo clinical supplement: {path}")
    clinical = pd.read_excel(path, sheet_name="clinical annotation", header=1)
    clinical = clinical.rename(columns={"Name of Sample": "sample", "Primary/ Recurrence": "presentation"})
    clinical["sample"] = clinical["sample"].astype(str)
    clinical["sample_base"] = clinical["sample"].str.replace(r"/R\d+$", "", regex=True)
    clinical["condition"] = clinical["presentation"].map(condition_from_text)
    clinical["PFS_years"] = pd.to_numeric(clinical["PFS [Years]"], errors="coerce")
    clinical["OS_years"] = pd.to_numeric(clinical["OS [years]"], errors="coerce")
    keep = ["sample_base", "Molecular group", "condition", "presentation", "PFS_years", "OS_years", "outcome"]
    collapsed = (
        clinical[keep]
        .sort_values("condition")
        .groupby("sample_base", as_index=False)
        .agg(
            molecular_group=("Molecular group", lambda values: next((v for v in values if pd.notna(v)), np.nan)),
            condition=("condition", lambda values: "Recurrent" if (values == "Recurrent").any() else ("Primary" if (values == "Primary").any() else "Unknown")),
            presentation=("presentation", lambda values: "; ".join(sorted(set(str(v) for v in values if pd.notna(v))))),
            PFS_years=("PFS_years", "max"),
            OS_years=("OS_years", "max"),
            outcome=("outcome", lambda values: next((v for v in values if pd.notna(v)), np.nan)),
        )
    )
    return collapsed


def parse_lr_genes_from_existing_outputs() -> list[str]:
    genes = set(KEY_LR_GENES)
    stats_path = OUTPUT_DIR / "fig4_survival" / "GSE64415_bulk_key_LR_stats.csv"
    if stats_path.exists():
        stats = pd.read_csv(stats_path)
        genes.update(stats["gene"].astype(str))
    for csv_path in [
        OUTPUT_DIR / "fig5_readable" / "Fig5B_mal_to_micro_top_interactions.csv",
        OUTPUT_DIR / "fig5_readable" / "Fig5B_micro_to_mal_top_interactions.csv",
    ]:
        if not csv_path.exists():
            continue
        interactions = pd.read_csv(csv_path)
        for value in interactions["interaction_name"].astype(str):
            for token in re.split(r"[_\-]", value):
                if re.fullmatch(r"[A-Z][A-Z0-9]{1,}", token):
                    genes.add(token)
    excluded = {"PGE2", "GLU", "CHOL", "CHOLESTEROL"}
    return sorted(genes - excluded)


def read_gene_list(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def list_sample_files() -> list[tuple[Path, Path]]:
    meta_files = sorted(path for path in ROOT_DIR.rglob("*_meta.txt.gz") if not path.name.startswith("._"))
    pairs = []
    for meta_path in meta_files:
        cm_path = Path(str(meta_path).replace("_meta.txt.gz", "_cm.txt.gz"))
        if cm_path.exists():
            pairs.append((meta_path, cm_path))
    return pairs


def compartment(annotation: str, malignant: str) -> str:
    if malignant == "Malignant":
        return "Neoplastic"
    text = str(annotation).lower()
    if "microglia" in text or "myeloid" in text or "macroph" in text:
        return "Myeloid"
    if "t-cell" in text or "b-cell" in text or "lymph" in text:
        return "Lymphocytes"
    return "Other non-malignant"


def sample_key(path: Path) -> str:
    return path.name.replace("_200113lj_meta.txt.gz", "").replace("_200226lj_meta.txt.gz", "")


def stream_one_sample(
    meta_path: Path,
    cm_path: Path,
    sample_info: pd.Series,
    undiff_genes: set[str],
    diff_genes: set[str],
    lr_genes: set[str],
) -> tuple[pd.DataFrame, dict[str, float], dict[str, dict[str, float]]]:
    meta = pd.read_csv(meta_path, sep="\t")
    with gzip.open(cm_path, "rt") as handle:
        cell_ids = handle.readline().rstrip("\n").split("\t")
        if len(cell_ids) != len(meta):
            raise ValueError(f"{cm_path} has {len(cell_ids)} cells but metadata has {len(meta)} rows")
        meta = meta.copy()
        meta["cell_id"] = cell_ids
        meta["sample_file"] = sample_key(meta_path)
        meta["condition"] = sample_info["condition"]
        meta["condition_source"] = sample_info["condition_source"]
        meta["presentation"] = sample_info.get("presentation", np.nan)
        meta["molecular_group"] = sample_info.get("molecular_group", np.nan)
        meta["PFS_years"] = sample_info.get("PFS_years", np.nan)
        meta["OS_years"] = sample_info.get("OS_years", np.nan)
        meta["outcome"] = sample_info.get("outcome", np.nan)
        meta["compartment"] = [compartment(a, m) for a, m in zip(meta["annotation"], meta["malignant"], strict=False)]

        undiff_sum = np.zeros(len(meta), dtype=np.float64)
        diff_sum = np.zeros(len(meta), dtype=np.float64)
        matched_undiff = 0
        matched_diff = 0
        lr_sums: dict[str, dict[str, float]] = {}
        pseudo_sum: dict[str, float] = {}

        malignant_mask = (meta["malignant"].to_numpy() == "Malignant")
        for line_number, line in enumerate(handle, start=2):
            gene, values = line.rstrip("\n").split("\t", 1)
            row = parse_tsv_floats(values, dtype=np.float32)
            if row.size != len(meta):
                raise ValueError(f"{cm_path} line {line_number} has {row.size} values")
            log_expr = np.log1p(row.astype(np.float64))
            if malignant_mask.any():
                pseudo_sum[gene] = float(log_expr[malignant_mask].mean())
            if gene in undiff_genes:
                undiff_sum += log_expr
                matched_undiff += 1
            if gene in diff_genes:
                diff_sum += log_expr
                matched_diff += 1
            if gene in lr_genes:
                lr_sums[gene] = {
                    "gene": gene,
                    **{
                        f"{label}__mean_log1p_tpm": float(log_expr[(meta["compartment"] == label).to_numpy()].mean())
                        if (meta["compartment"] == label).any() else np.nan
                        for label in ["Neoplastic", "Myeloid", "Lymphocytes", "Other non-malignant"]
                    },
                    **{
                        f"{label}__pct_expressing": float((row[(meta["compartment"] == label).to_numpy()] > 0).mean())
                        if (meta["compartment"] == label).any() else np.nan
                        for label in ["Neoplastic", "Myeloid", "Lymphocytes", "Other non-malignant"]
                    },
                }

    meta["undiff_module_score"] = undiff_sum / max(matched_undiff, 1)
    meta["diff_module_score"] = diff_sum / max(matched_diff, 1)
    meta["trajectory_score"] = meta["undiff_module_score"] - meta["diff_module_score"]
    meta["matched_undiff_genes"] = matched_undiff
    meta["matched_diff_genes"] = matched_diff
    return meta, pseudo_sum, lr_sums


def build_sample_info(meta_paths: list[Path]) -> pd.DataFrame:
    clinical = load_gojo_clinical()
    rows = []
    for meta_path in meta_paths:
        meta = pd.read_csv(meta_path, sep="\t", nrows=1)
        sample = str(meta["sample"].iloc[0])
        sample_file = sample_key(meta_path)
        sample_base = re.sub(r"Nuc\d*$", "", sample_file)
        sample_base = re.sub(r"Nuc$", "", sample_base)
        sample_base = re.sub(r"Dia.*$", "", sample_base)
        sample_base = re.sub(r"Rec.*$", "", sample_base)
        sample_base = sample_base if sample_base else sample
        clin = clinical[clinical["sample_base"] == sample_base]
        condition = condition_from_sample(sample_file)
        source = "sample_name"
        if condition is None and not clin.empty:
            condition = clin["condition"].iloc[0]
            source = "clinical_table"
        if condition is None:
            condition = "Unknown"
            source = "missing"
        rows.append(
            {
                "sample": sample,
                "sample_file": sample_file,
                "sample_base": sample_base,
                "condition": condition,
                "condition_source": source,
                "subtype": meta["subtype"].iloc[0],
                "molecular_group": clin["molecular_group"].iloc[0] if not clin.empty else np.nan,
                "presentation": clin["presentation"].iloc[0] if not clin.empty else np.nan,
                "PFS_years": clin["PFS_years"].iloc[0] if not clin.empty else np.nan,
                "OS_years": clin["OS_years"].iloc[0] if not clin.empty else np.nan,
                "outcome": clin["outcome"].iloc[0] if not clin.empty else np.nan,
            }
        )
    info = pd.DataFrame(rows)
    info.to_csv(v2_file("fig4_gojo_scrna/GSE141460_sample_metadata.csv"), index=False)
    return info


def savefig(path_stub: Path, width: float, height: float) -> None:
    plt.gcf().set_size_inches(width, height)
    plt.tight_layout()
    plt.savefig(path_stub.with_suffix(".pdf"))
    plt.savefig(path_stub.with_suffix(".png"), dpi=220)
    plt.close()


def plot_composition(cells: pd.DataFrame) -> None:
    comp = (
        cells.groupby(["sample_file", "condition", "subtype", "annotation"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    comp["fraction"] = comp["n_cells"] / comp.groupby("sample_file")["n_cells"].transform("sum")
    comp.to_csv(v2_file("fig4_gojo_scrna/GSE141460_celltype_composition_by_sample.csv"), index=False)
    top = comp.groupby("annotation")["n_cells"].sum().sort_values(ascending=False).head(14).index
    comp["annotation_plot"] = np.where(comp["annotation"].isin(top), comp["annotation"], "Other")
    plot_data = comp.groupby(["sample_file", "condition", "subtype", "annotation_plot"], dropna=False)["fraction"].sum().reset_index()
    order = (
        plot_data[["sample_file", "condition", "subtype"]]
        .drop_duplicates()
        .assign(condition_order=lambda data: data["condition"].map({"Primary": 0, "Recurrent": 1, "Unknown": 2}).fillna(3))
        .sort_values(["condition_order", "subtype", "sample_file"])
    )
    pivot = plot_data.pivot_table(index="sample_file", columns="annotation_plot", values="fraction", fill_value=0)
    pivot = pivot.loc[order["sample_file"]]
    colors = sns.color_palette("tab20", n_colors=len(pivot.columns))
    ax = pivot.plot(kind="bar", stacked=True, color=colors, linewidth=0, width=0.86)
    ax.set_title("Gojo GSE141460 scRNA cell-state composition by sample")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Cell fraction")
    ax.legend(title="Annotation", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.tick_params(axis="x", rotation=90)
    savefig(v2_file("fig4_gojo_scrna/Fig4B_GSE141460_cellstate_composition_by_sample_scrna"), 13.5, 6.2)


def compare_groups(summary: pd.DataFrame, subset: str | None = None) -> dict[str, object]:
    data = summary.copy()
    label = "All subtypes"
    if subset is not None:
        data = data[data["subtype"] == subset].copy()
        label = f"{subset} only"
    primary = data.loc[data["condition"] == "Primary", "mean_trajectory_score"]
    recurrent = data.loc[data["condition"] == "Recurrent", "mean_trajectory_score"]
    p_value = np.nan
    if len(primary) > 0 and len(recurrent) > 0:
        p_value = mannwhitneyu(primary, recurrent, alternative="two-sided").pvalue
    return {
        "comparison": f"{label}: recurrent vs primary",
        "primary_n": len(primary),
        "recurrent_n": len(recurrent),
        "primary_median": primary.median(),
        "recurrent_median": recurrent.median(),
        "mannwhitney_p": p_value,
    }


def plot_trajectory(cells: pd.DataFrame) -> pd.DataFrame:
    malignant = cells[(cells["malignant"] == "Malignant") & cells["condition"].isin(["Primary", "Recurrent"])].copy()
    malignant.to_csv(v2_file("fig4_gojo_scrna/GSE141460_malignant_cell_trajectory_scores.csv"), index=False)
    summary = (
        malignant.groupby(["sample_file", "sample", "condition", "subtype"], dropna=False)
        .agg(
            n_malignant_cells=("cell_id", "size"),
            mean_trajectory_score=("trajectory_score", "mean"),
            median_trajectory_score=("trajectory_score", "median"),
            mean_undiff_module_score=("undiff_module_score", "mean"),
            mean_diff_module_score=("diff_module_score", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(v2_file("fig4_gojo_scrna/GSE141460_sample_trajectory_summary.csv"), index=False)
    stats = pd.DataFrame([compare_groups(summary), compare_groups(summary, "PF-A"), compare_groups(summary, "ST-RELA")])
    stats.to_csv(v2_file("fig4_gojo_scrna/GSE141460_trajectory_condition_stats.csv"), index=False)

    sns.set_theme(style="whitegrid", context="paper")
    plt.figure()
    ax = sns.boxplot(data=summary, x="condition", y="mean_trajectory_score", order=["Primary", "Recurrent"], color="white", width=0.42, fliersize=0)
    sns.stripplot(data=summary, x="condition", y="mean_trajectory_score", order=["Primary", "Recurrent"], hue="subtype", palette="Set2", size=5, jitter=0.12, ax=ax)
    all_p = stats.loc[stats["comparison"] == "All subtypes: recurrent vs primary", "mannwhitney_p"].iloc[0]
    ax.text(0.02, 0.98, f"All subtypes p={all_p:.3g}", transform=ax.transAxes, ha="left", va="top", fontsize=9)
    ax.set_title("Gojo GSE141460 malignant trajectory score by presentation")
    ax.set_xlabel("")
    ax.set_ylabel("Sample mean trajectory score")
    ax.legend(title="Subtype", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    savefig(v2_file("fig4_gojo_scrna/Fig4_GSE141460_malignant_trajectory_score_primary_recurrent_scrna"), 6.8, 4.8)

    pf_a = summary[summary["subtype"] == "PF-A"].copy()
    if pf_a["condition"].nunique() == 2:
        plt.figure()
        ax = sns.boxplot(data=pf_a, x="condition", y="mean_trajectory_score", order=["Primary", "Recurrent"], color="white", width=0.42, fliersize=0)
        sns.stripplot(data=pf_a, x="condition", y="mean_trajectory_score", order=["Primary", "Recurrent"], color="#66C2A5", size=5, jitter=0.1, ax=ax)
        pf_p = stats.loc[stats["comparison"] == "PF-A only: recurrent vs primary", "mannwhitney_p"].iloc[0]
        ax.text(0.02, 0.98, f"PF-A p={pf_p:.3g}", transform=ax.transAxes, ha="left", va="top", fontsize=9)
        ax.set_title("Gojo GSE141460 PF-A trajectory score")
        ax.set_xlabel("")
        ax.set_ylabel("Sample mean trajectory score")
        savefig(v2_file("fig4_gojo_scrna/Fig4_GSE141460_PFA_trajectory_score_primary_recurrent_scrna"), 5.2, 4.5)
    return summary


def differential_expression(pseudobulk: pd.DataFrame, summary: pd.DataFrame, subset: str | None = None) -> pd.DataFrame:
    meta = summary[["sample_file", "condition", "subtype"]].drop_duplicates().copy()
    if subset is not None:
        meta = meta[meta["subtype"] == subset].copy()
    meta = meta[meta["condition"].isin(["Primary", "Recurrent"])]
    data = pseudobulk[meta["sample_file"]]
    primary_idx = np.where(meta["condition"].to_numpy() == "Primary")[0]
    recurrent_idx = np.where(meta["condition"].to_numpy() == "Recurrent")[0]
    rows = []
    for gene, values in zip(data.index, data.to_numpy(dtype=np.float64), strict=False):
        primary = values[primary_idx]
        recurrent = values[recurrent_idx]
        if len(primary) == 0 or len(recurrent) == 0:
            continue
        t_p = np.nan
        mw_p = np.nan
        if len(primary) >= 2 and len(recurrent) >= 2:
            t_p = ttest_ind(recurrent, primary, equal_var=False, nan_policy="omit").pvalue
            mw_p = mannwhitneyu(recurrent, primary, alternative="two-sided").pvalue
        rows.append(
            {
                "gene": gene,
                "primary_n": len(primary),
                "recurrent_n": len(recurrent),
                "primary_mean_log1p_tpm": float(primary.mean()),
                "recurrent_mean_log1p_tpm": float(recurrent.mean()),
                "recurrent_minus_primary": float(recurrent.mean() - primary.mean()),
                "ttest_p": float(t_p) if not np.isnan(t_p) else np.nan,
                "mannwhitney_p": float(mw_p) if not np.isnan(mw_p) else np.nan,
            }
        )
    deg = pd.DataFrame(rows)
    if not deg.empty:
        deg["ttest_fdr"] = bh_adjust(deg["ttest_p"].to_numpy(dtype=np.float64))
        deg = deg.sort_values(["ttest_p", "recurrent_minus_primary"], ascending=[True, False])
    return deg


def plot_lr(lr_records: list[dict[str, object]]) -> None:
    lr = pd.DataFrame(lr_records)
    lr.to_csv(v2_file("fig4_gojo_scrna/GSE141460_key_LR_expression_by_condition_compartment.csv"), index=False)
    lr = lr[lr["condition"].isin(["Primary", "Recurrent"])]
    plot = (
        lr.groupby(["gene", "condition", "compartment"], dropna=False)["mean_log1p_tpm"]
        .mean()
        .reset_index()
    )
    preferred = [gene for gene in KEY_LR_GENES if gene in plot["gene"].unique()]
    plot = plot[plot["gene"].isin(preferred)].copy()
    plot["column"] = plot["condition"] + " " + plot["compartment"]
    column_order = [
        "Primary Neoplastic", "Recurrent Neoplastic",
        "Primary Myeloid", "Recurrent Myeloid",
        "Primary Lymphocytes", "Recurrent Lymphocytes",
        "Primary Other non-malignant", "Recurrent Other non-malignant",
    ]
    matrix = plot.pivot_table(index="gene", columns="column", values="mean_log1p_tpm", fill_value=0)
    matrix = matrix.reindex(index=preferred, columns=[col for col in column_order if col in matrix.columns])
    z_matrix = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    plt.figure(figsize=(9.2, 7.2))
    ax = sns.heatmap(z_matrix, cmap="vlag", center=0, linewidths=0.35, linecolor="white", cbar_kws={"label": "row z-score"})
    ax.set_title("Gojo GSE141460 key ligand/receptor expression by condition")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.xticks(rotation=35, ha="right")
    savefig(v2_file("fig4_gojo_scrna/Fig5DE_GSE141460_key_LR_expression_condition_compartment_scrna"), 9.2, 7.2)


def main() -> None:
    if not ROOT_DIR.exists():
        raise FileNotFoundError(f"Missing extracted Gojo TPM/meta directory: {ROOT_DIR}")
    # OUT_DIR is already created at import time via v2_dir().
    pairs = list_sample_files()
    sample_info = build_sample_info([meta for meta, _ in pairs])
    info_by_file = sample_info.set_index("sample_file")
    undiff_genes = set(read_gene_list(dep_file("fig3", "undiff_genes.txt")))
    diff_genes = set(read_gene_list(dep_file("fig3", "diff_genes.txt")))
    lr_genes = set(parse_lr_genes_from_existing_outputs())

    all_cells = []
    pseudobulk_records: dict[str, dict[str, float]] = {}
    lr_records = []
    for meta_path, cm_path in pairs:
        key = sample_key(meta_path)
        cells, pseudo, lr_sums = stream_one_sample(
            meta_path, cm_path, info_by_file.loc[key], undiff_genes, diff_genes, lr_genes
        )
        all_cells.append(cells)
        pseudobulk_records[key] = pseudo
        for gene, row in lr_sums.items():
            for compartment_name in ["Neoplastic", "Myeloid", "Lymphocytes", "Other non-malignant"]:
                lr_records.append(
                    {
                        "sample_file": key,
                        "sample": cells["sample"].iloc[0],
                        "condition": cells["condition"].iloc[0],
                        "subtype": cells["subtype"].iloc[0],
                        "gene": gene,
                        "compartment": compartment_name,
                        "mean_log1p_tpm": row.get(f"{compartment_name}__mean_log1p_tpm", np.nan),
                        "pct_expressing": row.get(f"{compartment_name}__pct_expressing", np.nan),
                    }
                )

    cells = pd.concat(all_cells, ignore_index=True)
    cells.to_csv(v2_file("fig4_gojo_scrna/GSE141460_cell_metadata_with_trajectory_scores.csv"), index=False)
    plot_composition(cells)
    summary = plot_trajectory(cells)

    pseudobulk = pd.DataFrame(pseudobulk_records).sort_index()
    pseudobulk.to_csv(v2_file("fig4_gojo_scrna/GSE141460_malignant_pseudobulk_mean_log1p_tpm.csv"))
    all_deg = differential_expression(pseudobulk, summary)
    all_deg.to_csv(v2_file("fig4_gojo_scrna/GSE141460_malignant_pseudobulk_DEG_all_subtypes.csv"), index=False)
    pfa_deg = differential_expression(pseudobulk, summary, "PF-A")
    pfa_deg.to_csv(v2_file("fig4_gojo_scrna/GSE141460_malignant_pseudobulk_DEG_PFA_only.csv"), index=False)
    recurrent_high = pfa_deg[pfa_deg["recurrent_minus_primary"] > 0].copy()
    if recurrent_high.empty:
        recurrent_high = all_deg[all_deg["recurrent_minus_primary"] > 0].copy()
    recurrent_high.head(200).to_csv(v2_file("fig4_gojo_scrna/GSE141460_recurrent_high_top200.csv"), index=False)
    v2_file("fig4_gojo_scrna/GSE141460_recurrent_high_top200_genes.txt").write_text(
        "\n".join(recurrent_high.head(200)["gene"].astype(str)) + "\n"
    )
    plot_lr(lr_records)
    print(
        "Gojo GSE141460 validation complete: "
        f"{cells['sample_file'].nunique()} sample files, {len(cells)} cells, "
        f"conditions={cells['condition'].value_counts().to_dict()}"
    )
    print(f"Outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()