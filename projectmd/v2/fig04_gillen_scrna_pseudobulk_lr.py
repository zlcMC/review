#!/usr/bin/env python3
"""GSE125969 pseudobulk DEG and key ligand/receptor expression panels."""

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
    KEY_LR_GENES_CORE,
    OUTPUT_DIR,
    bh_adjust,
    dep_file,
    parse_tsv_floats,
    v2_dir,
    v2_file,
)


CACHE_DIR = OUTPUT_DIR / "external_cache" / "GSE125969"   # large cache stays at v1 location
OUT_DIR = v2_dir("fig4_scrna")
COUNT_PATH = CACHE_DIR / "GSE125969_count_matrix.tsv.gz"
CELL_SCORE_PATH = dep_file("fig4_scrna", "GSE125969_cell_metadata_with_trajectory_scores.csv")

KEY_LR_GENES = list(KEY_LR_GENES_CORE)


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


def load_cells() -> pd.DataFrame:
    if not CELL_SCORE_PATH.exists():
        raise FileNotFoundError(
            f"Missing {CELL_SCORE_PATH}. Run projectmd/v2/fig04_gillen_scrna_validation.py "
            "(or the legacy projectmd/fig04_gillen_scrna_validation.py) first."
        )
    cells = pd.read_csv(CELL_SCORE_PATH)
    cells["patient_id"] = cells["patient_id"].astype(int)
    cells["condition"] = cells["condition"].astype(str)
    cells["is_neoplastic"] = cells["is_neoplastic"].astype(bool)
    return cells


def make_group_codes(meta_ordered: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    neoplastic = meta_ordered["is_neoplastic"] & meta_ordered["condition"].isin(["Primary", "Recurrent"])
    group_table = (
        meta_ordered.loc[neoplastic, ["patient_id", "condition", "tumor_subtype"]]
        .drop_duplicates()
        .sort_values(["condition", "tumor_subtype", "patient_id"])
        .reset_index(drop=True)
    )
    group_table["group_id"] = np.arange(len(group_table))
    key_to_id = {
        (row.patient_id, row.condition, row.tumor_subtype): int(row.group_id)
        for row in group_table.itertuples(index=False)
    }
    codes = np.full(len(meta_ordered), -1, dtype=np.int32)
    for index, row in meta_ordered.loc[neoplastic, ["patient_id", "condition", "tumor_subtype"]].iterrows():
        codes[index] = key_to_id[(row["patient_id"], row["condition"], row["tumor_subtype"])]
    return codes, group_table


def stream_counts(cells: pd.DataFrame, lr_genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray], np.ndarray]:
    if not COUNT_PATH.exists():
        raise FileNotFoundError(f"Missing count matrix: {COUNT_PATH}")

    lr_gene_set = set(lr_genes)
    lr_counts: dict[str, np.ndarray] = {}
    genes: list[str] = []
    pseudobulk_rows: list[np.ndarray] = []

    with gzip.open(COUNT_PATH, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        meta_ordered = cells.set_index("cell_id").loc[cell_ids].reset_index()
        codes, group_table = make_group_codes(meta_ordered)
        neoplastic_mask = codes >= 0
        totals = np.zeros(len(cell_ids), dtype=np.float64)

        for line_number, line in enumerate(handle, start=2):
            gene, values = line.rstrip("\n").split("\t", 1)
            row = parse_tsv_floats(values, dtype=np.float32)
            if row.size != len(cell_ids):
                raise ValueError(f"Line {line_number} has {row.size} values, expected {len(cell_ids)}")
            totals += row

            pseudobulk = np.bincount(codes[neoplastic_mask], weights=row[neoplastic_mask], minlength=len(group_table))
            if pseudobulk.sum() > 0:
                genes.append(gene)
                pseudobulk_rows.append(pseudobulk.astype(np.float32))
            if gene in lr_gene_set:
                lr_counts[gene] = row.copy()

    pseudobulk_counts = pd.DataFrame(
        np.vstack(pseudobulk_rows),
        index=genes,
        columns=[f"{row.patient_id}|{row.condition}|{row.tumor_subtype}" for row in group_table.itertuples(index=False)],
    )
    group_table["sample"] = pseudobulk_counts.columns
    group_table.to_csv(v2_file("fig4_scrna/GSE125969_pseudobulk_sample_metadata.csv"), index=False)
    pseudobulk_counts.to_csv(v2_file("fig4_scrna/GSE125969_neoplastic_pseudobulk_counts.csv"))
    return pseudobulk_counts, group_table, lr_counts, totals


def differential_expression(counts: pd.DataFrame, groups: pd.DataFrame, subset: str | None = None) -> pd.DataFrame:
    selected_groups = groups.copy()
    if subset is not None:
        selected_groups = selected_groups[selected_groups["tumor_subtype"] == subset].copy()
    selected_counts = counts[selected_groups["sample"]]
    lib_sizes = selected_counts.sum(axis=0).to_numpy(dtype=np.float64)
    log_cpm = np.log2((selected_counts.to_numpy(dtype=np.float64) + 0.5) / (lib_sizes + 1.0) * 1e6)
    primary_idx = np.where(selected_groups["condition"].to_numpy() == "Primary")[0]
    recurrent_idx = np.where(selected_groups["condition"].to_numpy() == "Recurrent")[0]

    rows: list[dict[str, float | str | int]] = []
    for gene, values in zip(selected_counts.index, log_cpm, strict=False):
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
                "primary_mean_logCPM": float(primary.mean()),
                "recurrent_mean_logCPM": float(recurrent.mean()),
                "recurrent_minus_primary_logFC": float(recurrent.mean() - primary.mean()),
                "ttest_p": float(t_p) if not np.isnan(t_p) else np.nan,
                "mannwhitney_p": float(mw_p) if not np.isnan(mw_p) else np.nan,
            }
        )

    deg = pd.DataFrame(rows)
    if not deg.empty:
        deg["ttest_fdr"] = bh_adjust(deg["ttest_p"].to_numpy(dtype=np.float64))
        deg = deg.sort_values(["ttest_p", "recurrent_minus_primary_logFC"], ascending=[True, False])
    return deg


def save_deg_outputs(counts: pd.DataFrame, groups: pd.DataFrame) -> None:
    all_deg = differential_expression(counts, groups)
    all_deg.to_csv(v2_file("fig4_scrna/GSE125969_neoplastic_pseudobulk_DEG_all_subtypes.csv"), index=False)

    st_rela_deg = differential_expression(counts, groups, subset="ST-RELA")
    st_rela_deg.to_csv(v2_file("fig4_scrna/GSE125969_neoplastic_pseudobulk_DEG_ST_RELA_only.csv"), index=False)

    recurrent_high = all_deg[all_deg["recurrent_minus_primary_logFC"] > 0].copy()
    recurrent_high.head(200).to_csv(
        v2_file("fig4_scrna/GSE125969_neoplastic_pseudobulk_recurrent_high_top200.csv"), index=False
    )
    v2_file("fig4_scrna/GSE125969_neoplastic_pseudobulk_recurrent_high_top200_genes.txt").write_text(
        "\n".join(recurrent_high.head(200)["gene"].astype(str)) + "\n"
    )


def save_lr_expression(cells: pd.DataFrame, lr_counts: dict[str, np.ndarray], totals: np.ndarray) -> None:
    if not lr_counts:
        raise ValueError("No key LR genes were found in GSE125969 count matrix")

    with gzip.open(COUNT_PATH, "rt") as handle:
        cell_ids = handle.readline().rstrip("\n").split("\t")[1:]
    meta_ordered = cells.set_index("cell_id").loc[cell_ids].reset_index()
    meta_ordered["compartment"] = np.where(
        meta_ordered["is_neoplastic"],
        "Neoplastic",
        np.where(meta_ordered["cell_type"].isin(["Myeloid", "Lymphocytes"]), meta_ordered["cell_type"], "Other non-malignant"),
    )
    keep = meta_ordered["condition"].isin(["Primary", "Recurrent"])
    meta_kept = meta_ordered.loc[keep].copy()
    totals_kept = np.maximum(totals[keep.to_numpy()], 1.0)

    records: list[dict[str, object]] = []
    for gene, raw in lr_counts.items():
        expr = np.log1p(raw[keep.to_numpy()] / totals_kept * 1e4)
        expressed = raw[keep.to_numpy()] > 0
        tmp = meta_kept[["condition", "compartment"]].copy()
        tmp["expr"] = expr
        tmp["expressed"] = expressed
        grouped = tmp.groupby(["condition", "compartment"], dropna=False)
        for (condition, compartment), group in grouped:
            records.append(
                {
                    "gene": gene,
                    "condition": condition,
                    "compartment": compartment,
                    "n_cells": len(group),
                    "mean_log1p_cpm": group["expr"].mean(),
                    "pct_expressing": group["expressed"].mean(),
                }
            )

    lr_expr = pd.DataFrame(records)
    lr_expr.to_csv(v2_file("fig4_scrna/GSE125969_key_LR_expression_by_condition_compartment.csv"), index=False)

    preferred_genes = [gene for gene in KEY_LR_GENES if gene in lr_expr["gene"].unique()]
    lr_plot = lr_expr[lr_expr["gene"].isin(preferred_genes)].copy()
    lr_plot["column"] = lr_plot["condition"] + " " + lr_plot["compartment"]
    column_order = [
        "Primary Neoplastic", "Recurrent Neoplastic",
        "Primary Myeloid", "Recurrent Myeloid",
        "Primary Lymphocytes", "Recurrent Lymphocytes",
        "Primary Other non-malignant", "Recurrent Other non-malignant",
    ]
    matrix = lr_plot.pivot_table(index="gene", columns="column", values="mean_log1p_cpm", fill_value=0)
    matrix = matrix.reindex(index=preferred_genes, columns=[col for col in column_order if col in matrix.columns])
    z_matrix = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, np.nan), axis=0).fillna(0)

    sns.set_theme(style="white", context="paper")
    plt.figure(figsize=(9.2, 7.2))
    ax = sns.heatmap(z_matrix, cmap="vlag", center=0, linewidths=0.35, linecolor="white", cbar_kws={"label": "row z-score"})
    ax.set_title("GSE125969 key ligand/receptor expression by condition")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(v2_file("fig4_scrna/Fig5DE_GSE125969_key_LR_expression_condition_compartment_scrna.pdf"))
    plt.savefig(v2_file("fig4_scrna/Fig5DE_GSE125969_key_LR_expression_condition_compartment_scrna.png"), dpi=220)
    plt.close()


def main() -> None:
    # OUT_DIR is already created at import time via v2_dir().
    cells = load_cells()
    lr_genes = parse_lr_genes_from_existing_outputs()
    counts, groups, lr_counts, totals = stream_counts(cells, lr_genes)
    save_deg_outputs(counts, groups)
    save_lr_expression(cells, lr_counts, totals)
    print(
        "GSE125969 pseudobulk/LR complete: "
        f"{counts.shape[0]} genes x {counts.shape[1]} patient pseudobulks; "
        f"{len(lr_counts)} LR genes found."
    )
    print(f"Outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()