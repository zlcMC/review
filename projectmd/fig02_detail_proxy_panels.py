#!/usr/bin/env python
"""Generate Fig. 2 detail proxy panels from available GTE009 data.

The paper's Fig. 2D-F requires gene-level inferCNV values and chromosome-level
CNV matrices. Those matrices are not present in the public h5ad/RDS files in
this workspace. This script therefore creates clearly labelled proxy panels
from the available expression matrix, cell-level CNV score, subclone labels,
CytoTRACE, and trajectory score.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.stats import mannwhitneyu, pearsonr
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workspace_paths import output_path, raw_data_path


FIG2 = Path(output_path("fig2"))
FIG2.mkdir(parents=True, exist_ok=True)

PALETTE_CELL = {
    "NSC-like": "#E64B35",
    "RGC-like": "#4DBBD5",
    "Ast-like": "#00A087",
    "Epe-like": "#3C5488",
    "OPC-like": "#F39B7F",
    "Oli-like": "#8491B4",
    "Neu-like": "#91D1C2",
    "OPC": "#B09C85",
    "Mic": "#7E6148",
}
PALETTE_SUB = {
    "Subclone_1": "#D55E00",
    "Subclone_2": "#0072B2",
}


def get_vector(adata: sc.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    values = adata[:, gene].X
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).reshape(-1)


def mean_expression(adata: sc.AnnData, genes: list[str]) -> np.ndarray:
    genes = [gene for gene in genes if gene in adata.var_names]
    if not genes:
        return np.full(adata.n_obs, np.nan)
    matrix = adata[:, genes].X
    if sparse.issparse(matrix):
        return np.asarray(matrix.mean(axis=1)).reshape(-1)
    return np.asarray(matrix).mean(axis=1)


def zscore_rows(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame.astype(float)
    return values.sub(values.mean(axis=1), axis=0).div(values.std(axis=1) + 1e-6, axis=0)


def save_overview(adata: sc.AnnData) -> None:
    coords = adata.obsm["X_umap"]
    obs = adata.obs.copy()
    obs["UMAP_1"] = coords[:, 0]
    obs["UMAP_2"] = coords[:, 1]
    obs["DYNC2H1"] = get_vector(adata, "DYNC2H1")

    fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    panels = [
        ("Final_cluster", "categorical", "Cell state"),
        ("CNV_Cluster", "categorical", "CNV subclone"),
        ("CytoTRACE", "continuous", "CytoTRACE"),
        ("trajectory_score", "continuous", "Trajectory score"),
    ]
    for ax, (column, kind, title) in zip(axes.flat, panels):
        data = obs.dropna(subset=[column]) if column == "CNV_Cluster" else obs
        if kind == "categorical":
            if column == "Final_cluster":
                palette = PALETTE_CELL
            else:
                palette = PALETTE_SUB
            for label, group in data.groupby(column, observed=True):
                ax.scatter(group["UMAP_1"], group["UMAP_2"], s=3, lw=0,
                           color=palette.get(str(label), "#999999"), label=str(label), alpha=0.85)
            ax.legend(markerscale=4, fontsize=6, frameon=False, loc="best")
        else:
            sca = ax.scatter(data["UMAP_1"], data["UMAP_2"], c=data[column], s=3, lw=0,
                             cmap="viridis", alpha=0.9)
            fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("Fig. 2A proxy: GTE009 malignant-state overview", fontsize=12)
    fig.savefig(FIG2 / "Fig2A_GTE009_overview_proxy.pdf")
    fig.savefig(FIG2 / "Fig2A_GTE009_overview_proxy.png", dpi=220)
    plt.close(fig)


def save_dync2h1_proxy(adata: sc.AnnData, deg: pd.DataFrame) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs["DYNC2H1"] = get_vector(adata, "DYNC2H1")
    obs = obs[obs["CNV_Cluster"].isin(PALETTE_SUB)].copy()

    stats_rows = []
    for variable in ["DYNC2H1", "CNV_level", "CytoTRACE", "trajectory_score"]:
        left = obs.loc[obs["CNV_Cluster"] == "Subclone_1", variable].dropna()
        right = obs.loc[obs["CNV_Cluster"] == "Subclone_2", variable].dropna()
        stat = mannwhitneyu(left, right, alternative="two-sided") if len(left) and len(right) else None
        stats_rows.append({
            "comparison": "Subclone_1_vs_Subclone_2",
            "variable": variable,
            "Subclone_1_mean": float(left.mean()) if len(left) else np.nan,
            "Subclone_2_mean": float(right.mean()) if len(right) else np.nan,
            "mannwhitney_p": float(stat.pvalue) if stat else np.nan,
        })

    deg = deg.copy()
    deg["minus_log10_padj"] = -np.log10(deg["p_val_adj"].clip(lower=1e-300))
    deg["highlight"] = np.where(deg["gene"].eq("DYNC2H1"), "DYNC2H1", "Other DEG")

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    sns.scatterplot(
        data=deg, x="avg_log2FC", y="minus_log10_padj", hue="highlight",
        palette={"DYNC2H1": "#D55E00", "Other DEG": "#BDBDBD"},
        s=14, linewidth=0, ax=axes[0, 0], legend=False,
    )
    dync = deg[deg["gene"] == "DYNC2H1"]
    if not dync.empty:
        axes[0, 0].text(float(dync["avg_log2FC"].iloc[0]), float(dync["minus_log10_padj"].iloc[0]),
                        " DYNC2H1", fontsize=8, va="center")
    axes[0, 0].axvline(0, color="#555555", lw=0.7)
    axes[0, 0].set_title("Expression DEG volcano (Sub1 vs Sub2)")
    axes[0, 0].set_xlabel("avg log2FC")
    axes[0, 0].set_ylabel("-log10 adjusted P")

    sns.violinplot(data=obs, x="CNV_Cluster", y="DYNC2H1", hue="CNV_Cluster",
                   palette=PALETTE_SUB, inner="box", cut=0, legend=False, ax=axes[0, 1])
    axes[0, 1].set_title("DYNC2H1 expression by subclone")
    axes[0, 1].set_xlabel("")
    axes[0, 1].set_ylabel("log-normalized expression")

    sns.violinplot(data=obs, x="CNV_Cluster", y="CNV_level", hue="CNV_Cluster",
                   palette=PALETTE_SUB, inner="box", cut=0, legend=False, ax=axes[1, 0])
    axes[1, 0].set_title("Cell-level CNV score by subclone")
    axes[1, 0].set_xlabel("")
    axes[1, 0].set_ylabel("CNV_level")

    sns.violinplot(data=obs, x="CNV_Cluster", y="CytoTRACE", hue="CNV_Cluster",
                   palette=PALETTE_SUB, inner="box", cut=0, legend=False, ax=axes[1, 1])
    axes[1, 1].set_title("CytoTRACE by subclone")
    axes[1, 1].set_xlabel("")
    axes[1, 1].set_ylabel("CytoTRACE")

    fig.suptitle("Fig. 2D-F proxy: expression and cell-level CNV summaries", fontsize=12)
    fig.savefig(FIG2 / "Fig2D_F_DYNC2H1_expression_CNV_proxy.pdf")
    fig.savefig(FIG2 / "Fig2D_F_DYNC2H1_expression_CNV_proxy.png", dpi=220)
    plt.close(fig)
    return pd.DataFrame(stats_rows)


def save_cilium_proxy(adata: sc.AnnData, go: pd.DataFrame, deg: pd.DataFrame) -> pd.DataFrame:
    cilium_terms = go[go["Description"].str.contains("cilium|microtubule", case=False, na=False)]
    cilium_genes = []
    for genes in cilium_terms["geneID"].dropna().head(5):
        cilium_genes.extend(str(genes).split("/"))
    cilium_genes = sorted(set(gene for gene in cilium_genes if gene in adata.var_names))

    deg_rank = deg.set_index("gene").reindex(cilium_genes).dropna(subset=["avg_log2FC"])
    selected = deg_rank.sort_values("avg_log2FC", ascending=False).head(35).index.tolist()
    if "DYNC2H1" in cilium_genes and "DYNC2H1" not in selected:
        selected = ["DYNC2H1"] + selected[:34]

    obs = adata.obs.copy()
    obs["cilium_module_score"] = mean_expression(adata, selected)
    obs = obs[obs["CNV_Cluster"].isin(PALETTE_SUB)].copy()
    obs["state_subclone"] = obs["Final_cluster"].astype(str) + " | " + obs["CNV_Cluster"].astype(str)

    expr = adata[:, selected].X
    if sparse.issparse(expr):
        expr = expr.toarray()
    expr = pd.DataFrame(expr, index=adata.obs_names, columns=selected).loc[obs.index]
    group_means = expr.groupby(obs["state_subclone"]).mean().T
    group_means = group_means.loc[:, sorted(group_means.columns)]

    fig, ax = plt.subplots(figsize=(max(8, group_means.shape[1] * 0.45), max(7, len(selected) * 0.22)))
    sns.heatmap(zscore_rows(group_means), cmap="RdBu_r", center=0, ax=ax,
                cbar_kws={"label": "row Z-score"})
    ax.set_title("Fig. 2G proxy: cilium-related gene expression")
    ax.set_xlabel("Cell state | subclone")
    ax.set_ylabel("Gene")
    fig.tight_layout()
    fig.savefig(FIG2 / "Fig2G_cilium_expression_heatmap_proxy.pdf")
    fig.savefig(FIG2 / "Fig2G_cilium_expression_heatmap_proxy.png", dpi=220)
    plt.close(fig)

    epc = obs[obs["Final_cluster"].isin(["Epe-like", "Ast-like", "RGC-like", "NSC-like"])].copy()
    corr_all = pearsonr(obs["cilium_module_score"], obs["CytoTRACE"])
    corr_epc = pearsonr(epc["cilium_module_score"], epc["CytoTRACE"])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    sns.scatterplot(data=obs, x="CytoTRACE", y="cilium_module_score", hue="CNV_Cluster",
                    palette=PALETTE_SUB, s=8, linewidth=0, alpha=0.65, ax=axes[0])
    axes[0].set_title(f"All malignant cells: r={corr_all.statistic:.2f}, p={corr_all.pvalue:.1e}")
    axes[0].set_xlabel("CytoTRACE")
    axes[0].set_ylabel("Cilium module expression")
    sns.scatterplot(data=epc, x="CytoTRACE", y="cilium_module_score", hue="Final_cluster",
                    palette=PALETTE_CELL, s=10, linewidth=0, alpha=0.7, ax=axes[1])
    axes[1].set_title(f"Selected malignant states: r={corr_epc.statistic:.2f}, p={corr_epc.pvalue:.1e}")
    axes[1].set_xlabel("CytoTRACE")
    axes[1].set_ylabel("Cilium module expression")
    fig.suptitle("Fig. 2H-I proxy: cilium expression vs undifferentiation score", fontsize=12)
    fig.savefig(FIG2 / "Fig2H_I_cilium_CytoTRACE_correlation_proxy.pdf")
    fig.savefig(FIG2 / "Fig2H_I_cilium_CytoTRACE_correlation_proxy.png", dpi=220)
    plt.close(fig)

    pd.DataFrame({"gene": selected}).to_csv(FIG2 / "Fig2_cilium_proxy_genes.csv", index=False)
    return pd.DataFrame([
        {"comparison": "cilium_module_vs_CytoTRACE", "subset": "all_malignant",
         "pearson_r": corr_all.statistic, "p_value": corr_all.pvalue, "n_cells": len(obs)},
        {"comparison": "cilium_module_vs_CytoTRACE", "subset": "selected_malignant_states",
         "pearson_r": corr_epc.statistic, "p_value": corr_epc.pvalue, "n_cells": len(epc)},
    ])


def main() -> None:
    sns.set_theme(style="white", context="paper")
    adata = sc.read_h5ad(raw_data_path("GTE009_final.h5ad"))
    deg = pd.read_csv(output_path("fig1_supp/GTE009_Sub1vs2_DEG.csv"))
    go = pd.read_csv(output_path("fig1_supp/Fig1I_GO_subclone1.csv"))

    save_overview(adata)
    dync_stats = save_dync2h1_proxy(adata, deg)
    corr_stats = save_cilium_proxy(adata, go, deg)
    pd.concat([dync_stats, corr_stats], ignore_index=True).to_csv(
        FIG2 / "Fig2_detail_proxy_stats.csv", index=False
    )
    print("Wrote Fig. 2 proxy panels to", FIG2)


if __name__ == "__main__":
    main()
