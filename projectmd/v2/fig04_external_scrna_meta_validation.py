#!/usr/bin/env python3
"""Low-memory cross-cohort validation across Gillen and Gojo external scRNA outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import binomtest, mannwhitneyu, t

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import (  # noqa: E402
    OUTPUT_DIR,
    bh_adjust,
    dep_file,
    subtype_group,
    v2_dir,
)


OUT_DIR = v2_dir("fig4_external_scrna_meta")


def load_trajectory() -> pd.DataFrame:
    gillen = pd.read_csv(dep_file("fig4_scrna", "GSE125969_patient_trajectory_summary.csv"))
    gillen = gillen.rename(
        columns={
            "patient_id": "sample_id",
            "tumor_subtype": "subtype",
            "n_neoplastic_cells": "n_cells",
        }
    )
    gillen["cohort"] = "Gillen_GSE125969"

    gojo = pd.read_csv(dep_file("fig4_gojo_scrna", "GSE141460_sample_trajectory_summary.csv"))
    gojo = gojo.rename(columns={"sample_file": "sample_id", "n_malignant_cells": "n_cells"})
    gojo["cohort"] = "Gojo_GSE141460"

    keep = [
        "cohort",
        "sample_id",
        "condition",
        "subtype",
        "n_cells",
        "mean_trajectory_score",
        "median_trajectory_score",
        "mean_undiff_module_score",
        "mean_diff_module_score",
    ]
    combined = pd.concat([gillen[keep], gojo[keep]], ignore_index=True)
    combined = combined[combined["condition"].isin(["Primary", "Recurrent"])].copy()
    combined["subtype_group"] = combined["subtype"].map(subtype_group)
    combined["condition_recurrent"] = (combined["condition"] == "Recurrent").astype(float)
    return combined


def linear_model_condition(data: pd.DataFrame, covariates: list[str]) -> dict[str, object]:
    model_data = data.dropna(subset=["median_trajectory_score", "condition_recurrent", *covariates]).copy()
    y = model_data["median_trajectory_score"].to_numpy(dtype=float)
    design_parts = [pd.Series(1.0, index=model_data.index, name="Intercept"), model_data["condition_recurrent"]]
    for covariate in covariates:
        dummies = pd.get_dummies(model_data[covariate], prefix=covariate, drop_first=True, dtype=float)
        design_parts.append(dummies)
    design = pd.concat(design_parts, axis=1)
    design = design.loc[:, design.nunique(dropna=False) > 1]
    if "Intercept" not in design.columns:
        design.insert(0, "Intercept", 1.0)
    if "condition_recurrent" not in design.columns or len(model_data) < 4:
        return {
            "n": len(model_data),
            "rank": np.nan,
            "df_resid": np.nan,
            "condition_beta": np.nan,
            "condition_p": np.nan,
            "note": "condition term unavailable",
        }
    x = design.to_numpy(dtype=float)
    rank = int(np.linalg.matrix_rank(x))
    df_resid = int(len(model_data) - rank)
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    condition_index = list(design.columns).index("condition_recurrent")
    if df_resid <= 0:
        return {
            "n": len(model_data),
            "rank": rank,
            "df_resid": df_resid,
            "condition_beta": float(beta[condition_index]),
            "condition_p": np.nan,
            "note": "no residual degrees of freedom",
        }
    residual = y - x @ beta
    mse = float(np.sum(residual**2) / df_resid)
    covariance = mse * np.linalg.pinv(x.T @ x)
    se = float(np.sqrt(max(covariance[condition_index, condition_index], 0)))
    if se == 0:
        p_value = np.nan
    else:
        statistic = float(beta[condition_index] / se)
        p_value = float(2 * t.sf(abs(statistic), df_resid))
    return {
        "n": len(model_data),
        "rank": rank,
        "df_resid": df_resid,
        "condition_beta": float(beta[condition_index]),
        "condition_p": p_value,
        "note": "+".join(["condition", *covariates]),
    }


def mann_whitney_summary(data: pd.DataFrame) -> dict[str, object]:
    primary = data.loc[data["condition"] == "Primary", "median_trajectory_score"].dropna()
    recurrent = data.loc[data["condition"] == "Recurrent", "median_trajectory_score"].dropna()
    p_value = np.nan
    if len(primary) > 0 and len(recurrent) > 0:
        p_value = float(mannwhitneyu(primary, recurrent, alternative="two-sided").pvalue)
    return {
        "primary_n": len(primary),
        "recurrent_n": len(recurrent),
        "primary_median": float(primary.median()) if len(primary) else np.nan,
        "recurrent_median": float(recurrent.median()) if len(recurrent) else np.nan,
        "mannwhitney_p": p_value,
    }


def trajectory_stats(trajectory: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    subsets = {
        "all_samples": trajectory,
        "PF-A_only": trajectory[trajectory["subtype_group"] == "PF-A"],
        "ST-RELA_only": trajectory[trajectory["subtype_group"] == "ST-RELA"],
    }
    for subset_name, subset in subsets.items():
        base = {"subset": subset_name, **mann_whitney_summary(subset)}
        rows.append({"model": "pooled_mannwhitney", **base})
        if subset["cohort"].nunique() > 1:
            lm = linear_model_condition(subset, ["cohort"])
            rows.append({"model": "linear_condition_plus_cohort", **base, **lm})
        if subset["cohort"].nunique() > 1 and subset["subtype_group"].nunique() > 1:
            lm = linear_model_condition(subset, ["cohort", "subtype_group"])
            rows.append({"model": "linear_condition_plus_cohort_subtype", **base, **lm})
    return pd.DataFrame(rows)


def save_trajectory_plot(trajectory: pd.DataFrame) -> None:
    plot_data = trajectory.copy()
    plot_data["condition"] = pd.Categorical(plot_data["condition"], ["Primary", "Recurrent"], ordered=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    sns.boxplot(
        data=plot_data,
        x="condition",
        y="median_trajectory_score",
        hue="cohort",
        width=0.6,
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=plot_data,
        x="condition",
        y="median_trajectory_score",
        hue="cohort",
        dodge=True,
        linewidth=0.5,
        edgecolor="black",
        alpha=0.85,
        size=4,
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles, strict=False))
    ax.legend(unique.values(), unique.keys(), frameon=False, loc="best", title="Cohort")
    ax.axhline(0, color="0.6", linewidth=0.8, linestyle="--")
    ax.set_xlabel("")
    ax.set_ylabel("Median trajectory score")
    ax.set_title("External scRNA trajectory validation")
    sns.despine(fig=fig)
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig4_external_scRNA_combined_trajectory_score.{suffix}", dpi=300)
    plt.close(fig)


def load_deg(path: Path, cohort: str) -> pd.DataFrame:
    deg = pd.read_csv(path)
    delta_col = "recurrent_minus_primary_logFC" if "recurrent_minus_primary_logFC" in deg.columns else "recurrent_minus_primary"
    mean_primary_col = next(col for col in deg.columns if col.startswith("primary_mean"))
    mean_recurrent_col = next(col for col in deg.columns if col.startswith("recurrent_mean"))
    return pd.DataFrame(
        {
            "gene": deg["gene"].astype(str),
            f"{cohort}_primary_mean": deg[mean_primary_col],
            f"{cohort}_recurrent_mean": deg[mean_recurrent_col],
            f"{cohort}_delta": deg[delta_col],
            f"{cohort}_mannwhitney_p": deg.get("mannwhitney_p", np.nan),
            f"{cohort}_ttest_fdr": deg.get("ttest_fdr", np.nan),
        }
    )


def deg_concordance() -> tuple[pd.DataFrame, pd.DataFrame]:
    gillen = load_deg(dep_file("fig4_scrna", "GSE125969_neoplastic_pseudobulk_DEG_all_subtypes.csv"), "Gillen")
    gojo = load_deg(dep_file("fig4_gojo_scrna", "GSE141460_malignant_pseudobulk_DEG_all_subtypes.csv"), "Gojo")
    merged = gillen.merge(gojo, on="gene", how="inner")
    merged["Gillen_sign"] = np.sign(merged["Gillen_delta"])
    merged["Gojo_sign"] = np.sign(merged["Gojo_delta"])
    nonzero = merged[(merged["Gillen_sign"] != 0) & (merged["Gojo_sign"] != 0)].copy()
    nonzero["direction_concordant"] = nonzero["Gillen_sign"] == nonzero["Gojo_sign"]
    nonzero["both_recurrent_high"] = (nonzero["Gillen_delta"] > 0) & (nonzero["Gojo_delta"] > 0)
    nonzero["Gillen_recurrent_high_rank"] = (-nonzero["Gillen_delta"]).rank(method="average")
    nonzero["Gojo_recurrent_high_rank"] = (-nonzero["Gojo_delta"]).rank(method="average")
    nonzero["mean_recurrent_high_rank"] = nonzero[["Gillen_recurrent_high_rank", "Gojo_recurrent_high_rank"]].mean(axis=1)
    both_high = nonzero[nonzero["both_recurrent_high"]].sort_values("mean_recurrent_high_rank")
    return nonzero.sort_values("gene"), both_high


def save_deg_plot(deg: pd.DataFrame, lr_genes: set[str]) -> None:
    plot_data = deg.replace([np.inf, -np.inf], np.nan).dropna(subset=["Gillen_delta", "Gojo_delta"]).copy()
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    ax.scatter(plot_data["Gillen_delta"], plot_data["Gojo_delta"], s=8, alpha=0.18, color="#3b4252", linewidths=0)
    highlighted = plot_data[plot_data["gene"].isin(lr_genes)]
    ax.scatter(highlighted["Gillen_delta"], highlighted["Gojo_delta"], s=22, alpha=0.9, color="#b23a48", linewidths=0)
    for _, row in highlighted.sort_values("mean_recurrent_high_rank").head(12).iterrows():
        ax.text(row["Gillen_delta"], row["Gojo_delta"], row["gene"], fontsize=7, alpha=0.9)
    ax.axhline(0, color="0.45", linewidth=0.8)
    ax.axvline(0, color="0.45", linewidth=0.8)
    ax.set_xlabel("Gillen recurrent - primary logFC")
    ax.set_ylabel("Gojo recurrent - primary log1p TPM")
    ax.set_title("External scRNA DEG direction concordance")
    sns.despine(fig=fig)
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig4D_external_scRNA_DEG_direction_concordance.{suffix}", dpi=300)
    plt.close(fig)


def load_lr_delta() -> pd.DataFrame:
    gillen = pd.read_csv(dep_file("fig4_scrna", "GSE125969_key_LR_expression_by_condition_compartment.csv"))
    gojo = pd.read_csv(dep_file("fig4_gojo_scrna", "GSE141460_key_LR_expression_by_condition_compartment.csv"))

    gillen_wide = gillen.pivot_table(
        index=["gene", "compartment"],
        columns="condition",
        values="mean_log1p_cpm",
        aggfunc="mean",
    ).reset_index()
    gillen_wide["Gillen_delta"] = gillen_wide.get("Recurrent", np.nan) - gillen_wide.get("Primary", np.nan)

    gojo_mean = (
        gojo.groupby(["gene", "compartment", "condition"], as_index=False)["mean_log1p_tpm"]
        .mean()
    )
    gojo_wide = gojo_mean.pivot_table(
        index=["gene", "compartment"],
        columns="condition",
        values="mean_log1p_tpm",
        aggfunc="mean",
    ).reset_index()
    gojo_wide["Gojo_delta"] = gojo_wide.get("Recurrent", np.nan) - gojo_wide.get("Primary", np.nan)

    merged = gillen_wide[["gene", "compartment", "Gillen_delta"]].merge(
        gojo_wide[["gene", "compartment", "Gojo_delta"]],
        on=["gene", "compartment"],
        how="inner",
    )
    merged["Gillen_sign"] = np.sign(merged["Gillen_delta"])
    merged["Gojo_sign"] = np.sign(merged["Gojo_delta"])
    merged["direction_concordant"] = merged["Gillen_sign"] == merged["Gojo_sign"]
    merged["both_recurrent_high"] = (merged["Gillen_delta"] > 0) & (merged["Gojo_delta"] > 0)
    return merged.sort_values(["compartment", "gene"])


def save_lr_direction_plot(lr_delta: pd.DataFrame) -> None:
    selected = lr_delta[lr_delta["compartment"].isin(["Neoplastic", "Myeloid", "Lymphocytes"])].copy()
    selected = selected[selected["direction_concordant"] | selected["both_recurrent_high"]].copy()
    selected["label"] = selected["gene"] + " | " + selected["compartment"]
    selected["rank_value"] = selected[["Gillen_delta", "Gojo_delta"]].mean(axis=1)
    selected = selected.sort_values("rank_value", ascending=False).head(45)
    if selected.empty:
        return
    matrix = selected.set_index("label")[["Gillen_sign", "Gojo_sign"]]
    fig_height = max(4.5, 0.18 * len(matrix) + 1.8)
    fig, ax = plt.subplots(figsize=(4.8, fig_height))
    sns.heatmap(
        matrix,
        cmap=sns.color_palette(["#4c78a8", "#f2f2f2", "#b23a48"], as_cmap=True),
        vmin=-1,
        vmax=1,
        center=0,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Direction sign"},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Key LR recurrent-primary direction")
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(OUT_DIR / f"Fig5DE_external_scRNA_LR_direction_concordance.{suffix}", dpi=300)
    plt.close(fig)


def concordance_summary(deg: pd.DataFrame, lr_delta: pd.DataFrame) -> pd.DataFrame:
    deg_successes = int(deg["direction_concordant"].sum())
    deg_total = int(deg["direction_concordant"].notna().sum())
    deg_p = float(binomtest(deg_successes, deg_total, 0.5, alternative="greater").pvalue) if deg_total else np.nan

    lr_clean = lr_delta[(lr_delta["Gillen_sign"] != 0) & (lr_delta["Gojo_sign"] != 0)].copy()
    lr_successes = int(lr_clean["direction_concordant"].sum())
    lr_total = int(lr_clean["direction_concordant"].notna().sum())
    lr_p = float(binomtest(lr_successes, lr_total, 0.5, alternative="greater").pvalue) if lr_total else np.nan

    neo = lr_clean[lr_clean["compartment"] == "Neoplastic"]
    neo_successes = int(neo["direction_concordant"].sum())
    neo_total = int(neo["direction_concordant"].notna().sum())
    neo_p = float(binomtest(neo_successes, neo_total, 0.5, alternative="greater").pvalue) if neo_total else np.nan

    return pd.DataFrame(
        [
            {
                "analysis": "all_gene_DEG_direction",
                "concordant": deg_successes,
                "total": deg_total,
                "fraction_concordant": deg_successes / deg_total if deg_total else np.nan,
                "binomial_p_greater_than_half": deg_p,
            },
            {
                "analysis": "key_LR_direction_all_compartments",
                "concordant": lr_successes,
                "total": lr_total,
                "fraction_concordant": lr_successes / lr_total if lr_total else np.nan,
                "binomial_p_greater_than_half": lr_p,
            },
            {
                "analysis": "key_LR_direction_neoplastic",
                "concordant": neo_successes,
                "total": neo_total,
                "fraction_concordant": neo_successes / neo_total if neo_total else np.nan,
                "binomial_p_greater_than_half": neo_p,
            },
        ]
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")

    trajectory = load_trajectory()
    trajectory.to_csv(OUT_DIR / "external_scRNA_combined_trajectory_samples.csv", index=False)
    stats = trajectory_stats(trajectory)
    stats.to_csv(OUT_DIR / "external_scRNA_trajectory_meta_stats.csv", index=False)
    save_trajectory_plot(trajectory)

    deg, both_high = deg_concordance()
    deg.to_csv(OUT_DIR / "external_scRNA_DEG_direction_concordance.csv", index=False)
    both_high.head(200).to_csv(OUT_DIR / "external_scRNA_concordant_recurrent_high_top200.csv", index=False)

    lr_delta = load_lr_delta()
    lr_delta.to_csv(OUT_DIR / "external_scRNA_key_LR_delta_concordance.csv", index=False)
    lr_genes = set(lr_delta["gene"].astype(str))
    save_deg_plot(deg, lr_genes)
    save_lr_direction_plot(lr_delta)

    concordance = concordance_summary(deg, lr_delta)
    concordance.to_csv(OUT_DIR / "external_scRNA_direction_concordance_summary.csv", index=False)

    print(f"Wrote combined external scRNA validation outputs to {OUT_DIR}")
    print(stats.to_string(index=False))
    print(concordance.to_string(index=False))


if __name__ == "__main__":
    main()