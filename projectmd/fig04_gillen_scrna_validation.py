#!/usr/bin/env python3
"""External Gillen GSE125969 scRNA validation for Fig. 4B/4D-style panels."""

from __future__ import annotations

import gzip
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

from workspace_paths import OUTPUT_DIR, output_path


CACHE_DIR = OUTPUT_DIR / "external_cache" / "GSE125969"
PMC_DIR = OUTPUT_DIR / "external_cache" / "pmc_supp"
OUT_DIR = OUTPUT_DIR / "fig4_scrna"

METADATA_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125969/suppl/GSE125969_cell_metadata.tsv.gz"
COUNT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125969/suppl/GSE125969_count_matrix.tsv.gz"

NON_TUMOR_CELL_TYPES = {"Myeloid", "Lymphocytes", "Oligodendrocytes", "Doublets"}


def download_if_needed(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    print(f"Downloading {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)


def read_gene_list(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def load_clinical() -> pd.DataFrame:
    clinical_path = PMC_DIR / "Gillen_supplement_2.xlsx"
    if not clinical_path.exists():
        raise FileNotFoundError(f"Missing Gillen supplement: {clinical_path}")
    clinical = pd.read_excel(clinical_path, sheet_name="Patient and sample details", header=1)
    presentation_col = "presentation (1= presentation, 2= 1st recurrence)"
    clinical = clinical.rename(
        columns={
            "patient ID": "patient_id",
            presentation_col: "presentation",
            "PFS (months)": "PFS_time_months",
            "recurrence event": "PFS_event",
            "OS (months)": "OS_time_months",
            "death event": "OS_event",
            "WHO grade": "who_grade",
            "age at Dx (years)": "age_years",
        }
    )
    keep = [
        "patient_id",
        "presentation",
        "who_grade",
        "age_years",
        "sex",
        "PFS_time_months",
        "PFS_event",
        "OS_time_months",
        "OS_event",
    ]
    clinical = clinical[keep].copy()
    clinical["patient_id"] = clinical["patient_id"].astype("Int64")
    clinical["presentation"] = clinical["presentation"].astype("Int64")
    clinical["condition"] = clinical["presentation"].map({1: "Primary", 2: "Recurrent"})
    return clinical


def load_metadata(metadata_path: Path, clinical: pd.DataFrame) -> pd.DataFrame:
    meta = pd.read_csv(metadata_path, sep="\t")
    meta["patient_id"] = meta["cell_id"].str.extract(r"foreman_(\d+)").astype(int)
    meta["is_neoplastic"] = ~meta["cell_type"].isin(NON_TUMOR_CELL_TYPES)
    merged = meta.merge(clinical, on="patient_id", how="left")
    merged["condition"] = merged["condition"].fillna("Unknown")
    return merged


def stream_signature_scores(
    count_path: Path,
    undiff_genes: list[str],
    diff_genes: list[str],
) -> pd.DataFrame:
    target_genes = set(undiff_genes) | set(diff_genes)
    counts_by_gene: dict[str, np.ndarray] = {}

    with gzip.open(count_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        totals = np.zeros(len(cell_ids), dtype=np.float64)

        for line_number, line in enumerate(handle, start=2):
            gene, values = line.rstrip("\n").split("\t", 1)
            row = np.fromstring(values, sep="\t", dtype=np.float32)
            if row.size != len(cell_ids):
                raise ValueError(f"Line {line_number} has {row.size} values, expected {len(cell_ids)}")
            totals += row
            if gene in target_genes:
                if gene in counts_by_gene:
                    counts_by_gene[gene] += row
                else:
                    counts_by_gene[gene] = row.copy()

    totals = np.maximum(totals, 1.0)

    def module_score(genes: list[str]) -> tuple[np.ndarray, list[str]]:
        matched = [gene for gene in genes if gene in counts_by_gene]
        if not matched:
            raise ValueError("No signature genes matched the GSE125969 count matrix")
        matrix = np.vstack([counts_by_gene[gene] for gene in matched]).astype(np.float64)
        normalized = np.log1p(matrix / totals * 1e4)
        return normalized.mean(axis=0), matched

    undiff_score, matched_undiff = module_score(undiff_genes)
    diff_score, matched_diff = module_score(diff_genes)
    return pd.DataFrame(
        {
            "cell_id": cell_ids,
            "undiff_module_score": undiff_score,
            "diff_module_score": diff_score,
            "trajectory_score": undiff_score - diff_score,
            "matched_undiff_genes": len(matched_undiff),
            "matched_diff_genes": len(matched_diff),
        }
    )


def savefig(path_stub: Path, width: float = 8, height: float = 5) -> None:
    plt.gcf().set_size_inches(width, height)
    plt.tight_layout()
    plt.savefig(path_stub.with_suffix(".pdf"))
    plt.savefig(path_stub.with_suffix(".png"), dpi=220)
    plt.close()


def plot_celltype_composition(cells: pd.DataFrame) -> None:
    composition = (
        cells.groupby(["patient_id", "condition", "tumor_subtype", "cell_type"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    totals = composition.groupby("patient_id")["n_cells"].transform("sum")
    composition["fraction"] = composition["n_cells"] / totals
    composition.to_csv(output_path("fig4_scrna/GSE125969_celltype_composition_by_patient.csv"), index=False)

    top_cell_types = composition.groupby("cell_type")["n_cells"].sum().sort_values(ascending=False).head(12).index
    composition["cell_type_plot"] = np.where(composition["cell_type"].isin(top_cell_types), composition["cell_type"], "Other")
    plot_data = (
        composition.groupby(["patient_id", "condition", "tumor_subtype", "cell_type_plot"], dropna=False)["fraction"]
        .sum()
        .reset_index()
    )
    order = (
        plot_data[["patient_id", "condition", "tumor_subtype"]]
        .drop_duplicates()
        .assign(condition_order=lambda data: data["condition"].map({"Primary": 0, "Recurrent": 1}).fillna(2))
        .sort_values(["condition_order", "tumor_subtype", "patient_id"])
    )
    order["patient_label"] = order["patient_id"].astype(str) + "\n" + order["condition"].str[0]
    plot_data = plot_data.merge(order[["patient_id", "patient_label"]], on="patient_id", how="left")
    pivot = plot_data.pivot_table(index="patient_label", columns="cell_type_plot", values="fraction", fill_value=0)
    pivot = pivot.loc[order["patient_label"]]

    colors = dict(zip(pivot.columns, sns.color_palette("tab20", n_colors=len(pivot.columns))))
    ax = pivot.plot(kind="bar", stacked=True, color=[colors[col] for col in pivot.columns], linewidth=0, width=0.86)
    ax.set_title("GSE125969 scRNA cell-type composition by patient")
    ax.set_xlabel("Patient (P=primary, R=recurrent)")
    ax.set_ylabel("Cell fraction")
    ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.tick_params(axis="x", labelrotation=90)
    savefig(output_path("fig4_scrna/Fig4B_GSE125969_celltype_composition_by_patient_scrna"), width=12.5, height=6.4)

    condition_mean = (
        composition.groupby(["condition", "cell_type"], dropna=False)["n_cells"]
        .sum()
        .reset_index()
    )
    condition_mean["fraction"] = condition_mean["n_cells"] / condition_mean.groupby("condition")["n_cells"].transform("sum")
    condition_mean.to_csv(output_path("fig4_scrna/GSE125969_celltype_composition_by_condition.csv"), index=False)


def plot_trajectory_scores(cells: pd.DataFrame) -> None:
    neoplastic = cells[cells["is_neoplastic"] & cells["condition"].isin(["Primary", "Recurrent"])].copy()
    neoplastic.to_csv(output_path("fig4_scrna/GSE125969_neoplastic_cell_trajectory_scores.csv"), index=False)

    patient_summary = (
        neoplastic.groupby(["patient_id", "condition", "tumor_subtype"], dropna=False)
        .agg(
            n_neoplastic_cells=("cell_id", "size"),
            mean_trajectory_score=("trajectory_score", "mean"),
            median_trajectory_score=("trajectory_score", "median"),
            mean_undiff_module_score=("undiff_module_score", "mean"),
            mean_diff_module_score=("diff_module_score", "mean"),
        )
        .reset_index()
    )
    patient_summary.to_csv(output_path("fig4_scrna/GSE125969_patient_trajectory_summary.csv"), index=False)

    def comparison_row(data: pd.DataFrame, label: str) -> dict[str, object]:
        primary = data.loc[data["condition"] == "Primary", "mean_trajectory_score"]
        recurrent = data.loc[data["condition"] == "Recurrent", "mean_trajectory_score"]
        p_value = np.nan
        if len(primary) > 0 and len(recurrent) > 0:
            p_value = mannwhitneyu(primary, recurrent, alternative="two-sided").pvalue
        return {
            "comparison": label,
            "primary_n": len(primary),
            "recurrent_n": len(recurrent),
            "primary_median": primary.median(),
            "recurrent_median": recurrent.median(),
            "mannwhitney_p": p_value,
            "matched_undiff_genes": int(neoplastic["matched_undiff_genes"].iloc[0]),
            "matched_diff_genes": int(neoplastic["matched_diff_genes"].iloc[0]),
        }

    summary_stats = pd.DataFrame(
        [
            comparison_row(patient_summary, "All subtypes: recurrent vs primary"),
            comparison_row(
                patient_summary[patient_summary["tumor_subtype"] == "ST-RELA"],
                "ST-RELA only: recurrent vs primary",
            ),
        ]
    )
    summary_stats.to_csv(output_path("fig4_scrna/GSE125969_trajectory_condition_stats.csv"), index=False)

    primary = patient_summary.loc[patient_summary["condition"] == "Primary", "mean_trajectory_score"]
    recurrent = patient_summary.loc[patient_summary["condition"] == "Recurrent", "mean_trajectory_score"]
    p_value = summary_stats.loc[summary_stats["comparison"] == "All subtypes: recurrent vs primary", "mannwhitney_p"].iloc[0]

    sns.set_theme(style="whitegrid", context="paper")
    plt.figure()
    ax = sns.boxplot(
        data=patient_summary,
        x="condition",
        y="mean_trajectory_score",
        order=["Primary", "Recurrent"],
        color="white",
        width=0.42,
        fliersize=0,
        linewidth=1.1,
    )
    sns.stripplot(
        data=patient_summary,
        x="condition",
        y="mean_trajectory_score",
        order=["Primary", "Recurrent"],
        hue="tumor_subtype",
        palette="Set2",
        size=5,
        jitter=0.12,
        ax=ax,
    )
    ax.set_title("GSE125969 external scRNA trajectory score by presentation")
    ax.set_xlabel("")
    ax.set_ylabel("Patient mean trajectory score\nmean(undiff genes) - mean(diff genes)")
    ax.text(
        0.02,
        0.98,
        f"Primary n={len(primary)}, recurrent n={len(recurrent)}; Mann-Whitney p={p_value:.3g}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    ax.legend(title="Subtype", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    savefig(output_path("fig4_scrna/Fig4_GSE125969_neoplastic_trajectory_score_primary_recurrent_scrna"), width=6.8, height=4.8)

    st_rela = patient_summary[patient_summary["tumor_subtype"] == "ST-RELA"].copy()
    if st_rela["condition"].nunique() == 2:
        st_p_value = summary_stats.loc[
            summary_stats["comparison"] == "ST-RELA only: recurrent vs primary", "mannwhitney_p"
        ].iloc[0]
        plt.figure()
        ax = sns.boxplot(
            data=st_rela,
            x="condition",
            y="mean_trajectory_score",
            order=["Primary", "Recurrent"],
            color="white",
            width=0.42,
            fliersize=0,
            linewidth=1.1,
        )
        sns.stripplot(
            data=st_rela,
            x="condition",
            y="mean_trajectory_score",
            order=["Primary", "Recurrent"],
            color="#8DA0CB",
            size=5,
            jitter=0.08,
            ax=ax,
        )
        ax.set_title("GSE125969 ST-RELA scRNA trajectory score")
        ax.set_xlabel("")
        ax.set_ylabel("Patient mean trajectory score")
        ax.text(
            0.02,
            0.98,
            f"Primary n={(st_rela['condition'] == 'Primary').sum()}, recurrent n={(st_rela['condition'] == 'Recurrent').sum()}; p={st_p_value:.3g}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
        )
        savefig(output_path("fig4_scrna/Fig4_GSE125969_ST_RELA_trajectory_score_primary_recurrent_scrna"), width=5.2, height=4.5)


def main() -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    metadata_path = CACHE_DIR / "GSE125969_cell_metadata.tsv.gz"
    count_path = CACHE_DIR / "GSE125969_count_matrix.tsv.gz"
    download_if_needed(METADATA_URL, metadata_path)
    download_if_needed(COUNT_URL, count_path)

    clinical = load_clinical()
    meta = load_metadata(metadata_path, clinical)
    undiff_genes = read_gene_list(output_path("fig3/undiff_genes.txt"))
    diff_genes = read_gene_list(output_path("fig3/diff_genes.txt"))
    scores = stream_signature_scores(count_path, undiff_genes, diff_genes)
    cells = meta.merge(scores, on="cell_id", how="left")
    cells.to_csv(output_path("fig4_scrna/GSE125969_cell_metadata_with_trajectory_scores.csv"), index=False)

    plot_celltype_composition(cells)
    plot_trajectory_scores(cells)

    print(
        "GSE125969 validation complete: "
        f"{cells['patient_id'].nunique()} patients, {len(cells)} cells, "
        f"matched undiff={int(scores['matched_undiff_genes'].iloc[0])}, "
        f"diff={int(scores['matched_diff_genes'].iloc[0])}."
    )
    print(f"Outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()