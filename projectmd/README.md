# 复现 fimmu-13-903246

> Wu et al., *Frontiers in Immunology* 13:903246 (2022)
> "Immune Crosstalk in Relapsed Ependymoma" (室管膜瘤 scRNA-seq)
> 原文 PDF: [read/fimmu-13-903246.pdf](../read/fimmu-13-903246.pdf)

## 1. 数据

原始数据位于 `projectfile/`（只读）：

| 文件 | 说明 |
|---|---|
| `GTE001.h5ad` / `GTE002.h5ad` / `GTE009.h5ad` / `GTE012.h5ad` | 4 例 EPN 样本（与论文 Figure 1 对应） |
| `GTE009_clustered.h5ad` | GTE009 已聚类版本 |
| `GTE009_final.h5ad` | GTE009 最终注释版本 |
| `GTE009_seurat.rds` | GTE009 Seurat 对象 |

所有运行产物写到 `output/`，不要写回 `projectfile/` 或 `read/`。

## 2. 软件栈与环境

文章工作流以 **R/Seurat v3.2.0** 为主，部分分析在 Python。下面列两个分工。

### 2.1 Python 环境 `epn2`

- 位置：`/home/zlcmc/miniconda3/envs/epn2`（Win 路径 `\\wsl$\Ubuntu-24.04\home\zlcmc\miniconda3\envs\epn2`）
- Python 3.10.20
- Jupyter kernel 名称：`Python (epn2)`
- 重建脚本：[setup_install_python_extras.sh](setup_install_python_extras.sh)

| 包 | 版本 | 论文中用途 |
|---|---|---|
| scanpy | 1.11.5 | 通用 scRNA-seq 分析 / 数据读写 |
| anndata | 0.11.4 | h5ad I/O |
| numpy | <2 | velocyto 编译依赖 |
| leidenalg / python-igraph | 0.11 | 聚类 |
| harmonypy | 0.2.0 | 整合（备选） |
| scrublet | 0.2.3 | 双细胞检测（备选 DoubletFinder） |
| scanorama | 1.7.4 | 整合（备选） |
| matplotlib / seaborn | — | 绘图 |
| **scVelo** | 0.3.4 | RNA velocity（Figure 2B, Supp. Fig 6A） |
| **velocyto** | 0.17.17 | 由 BAM → spliced/unspliced loom |
| **pySCENIC** | 0.12.1 | 调控网络（Figure 6；Supp. Fig 9） |
| jupyter / ipykernel | — | Notebook |

### 2.2 R 环境 `epn2_r`

- 位置：`/home/zlcmc/miniconda3/envs/epn2_r`
- **R 4.4.3**（原文 3.6.3；因 Bioconductor 最新版要求较新 R）
- **Seurat 5.4.0**（原文 3.2.0；API 大致兼容，某些老包需 `options(Seurat.object.assay.version="v3")`）
- yml: [epn2_r.environment.yml](epn2_r.environment.yml)
- Jupyter kernel 名称：`R (epn2_r)`

已安装关键包：monocle 2.34.0, infercnv 1.22.0, slingshot 2.14.0, TSCAN 1.44.0, clusterProfiler 4.14.0, NMF 0.21, ComplexHeatmap 2.22, circlize 0.4.18, survminer 0.5.2, corrplot 0.95, IRkernel 1.3.2。

从 GitHub 安装（yml 无法覆盖）：

```bash
conda run -n epn2_r R -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'
conda run -n epn2_r R -e 'remotes::install_github("jinworks/CellChat")'
```

手动 / 替代安装：

- **CytoTRACE2**（官方升级版，替代原版 CytoTRACE）— 已装 `CytoTRACE2 1.1.0`（`remotes::install_github("digitalcytometry/cytotrace2", subdir="cytotrace2_r")`）
- **presto** — 未装：源码 C++11 与 RcppArmadillo 要求的 C++14 冲突，失败。CellChat 在没有 presto 时会自动 fallback 到 Seurat 自己的 DEG，不影响结果，只是慢一点
- **SCUBI** — 可视化，非必需，用到再装
- **CSOMAP** — 3D ligand-receptor，非必需，用到再装

## 2.3 数据状况（**重要**）

作者上传到 GEO 的 metadata 已经**包含完整的论文中间结果**：

| metadata 列 | 论文方法步骤 |
|---|---|
| `DF_hi.lo` | DoubletFinder 双细胞标记 |
| `S.Score` / `G2M.Score` / `Phase` | 细胞周期 |
| `CNV_level`, `CNV_level_ave`, `CNV_cluster`, `CNV_Cluster` | inferCNV + 亚克隆（Subclone_1/2/3） |
| `CytoTRACE` | 未分化打分 |
| `Step1_cell_type`, `Final_cluster`, `Brief_cluster` | 细胞类型注释（Signature Enrichment 结果） |

**结论**：inferCNV / CytoTRACE / DoubletFinder / 细胞类型注释 **不需要从头跑**，直接用 metadata 即可复现论文图表。从头重跑仅用于验证或灵敏度分析。

**无法跑的步骤**：

| 图 | 缺失数据 | 说明 |
|---|---|---|
| Figure 2B / Supp. Fig 6A — RNA velocity | 10x BAM 文件 | 需从 GEO 下 fastq → CellRanger → velocyto loom，约 100 GB / 数小时 |
| Supp. Fig 2 — WES CNVkit | WES fastq | WES 原始数据不在 projectfile/ |
| Figure 4 — 生存分析 | 外部 expression + OS/PFS | Gillen Table S1 已接入 GSE125861 primary 样本并输出 OS/PFS KM；外部 scRNA 验证已接入 Gillen GSE125969 与 Gojo GSE141460，结论按 supportive/proxy 表述 |

## 2.4 命令行工具（从 fastq 重处理才需要）

当前 `projectfile/` 已是 h5ad + 带注释的 csv，下列通常不需要：

---

文章用 R v3.6.3 + Seurat v3.2.0。以下包/用途：

| R 包 | 论文中用途 |
|---|---|
| Seurat 3.2.0 | 预处理、聚类、DEG（贯穿） |
| DoubletFinder | 双细胞剔除 |
| inferCNV | CNV 推断（Figure 1A, C, D） |
| CytoTRACE | 未分化打分（Figure 1B） |
| CellChat | 细胞通讯（Figure 4A-B） |
| Monocle (v2, DDRTree) | 原文拟时序方法（Figure 2C, Supp. Fig 6B）；当前 R 4.4 + dplyr/igraph 新环境不稳定，最终复现不再保留 Monocle2 运行脚本 |
| TSCAN | 拟时序 |
| **Slingshot** | 当前 Figure 2C 最终替代拟时序方案 |
| clusterProfiler | KEGG/GO 富集 |
| corrplot 0.90 | 相关性矩阵 |
| survminer / survival | 生存分析 |
| coin | 置换检验 |
| ggplot2 | 绘图 |
| SCUBI | 单细胞可视化 |
| CSOMAP | 3D 空间 ligand-receptor |

> 现有 `epn` 环境含 R（r-mutex）但未确认 Seurat 3.2.0；如需复用先 `conda run -n epn R -e 'packageVersion("Seurat")'` 检查。

### 2.4.1 命令行/外部工具（仅用于从 fastq 重处理）

> 当前 `projectfile/` 已是 h5ad，下列工具一般不需要执行。

- CellRanger v3.1.0 + 参考基因组 hg19-3.0.0
- BWA + Picard（WES 比对去重）
- CNVkit（WES CNV）
- GraphPad Prism v8.3.0（统计；可用 scipy/statsmodels 替代）

## 3. 复现计划（利用作者 GEO metadata）

如果要按依赖顺序复习或重跑，请先看 [RUN_ORDER.md](RUN_ORDER.md)。下面的表只说明每个脚本对应的论文图和当前状态。

| 脚本 | 环境 | 论文图表 | 状态 |
|---|---|---|---|
| `fig01_tsne_cnv_cytotrace_seurat.R` | epn2_r | Figure 1A-E（CNV/CytoTRACE/subclone/Brief_cluster on tSNE） | ✅ 已完成四个样本，输出 `output/fig1/GTE00*_tsne_panels.pdf` 与 Seurat RDS |
| `fig03_subclone_undiff_diff_score.R` | epn2_r | Figure 3（GTE009 subclone DEG → Undiff-Diff score） | ✅ 已写完，依赖 Figure 1 的 GTE009 Seurat 输出 |
| `fig02c_slingshot_final_panel.R` | epn2_r | Figure 2C（GTE009 恶性细胞拟时序；Slingshot 替代 Monocle2） | ✅ 已完成，主轨迹 NSC-like → RGC-like → Ast-like → Epe-like |
| `fig02_slingshot_basic.R` | epn2_r | Figure 2C 基础版 Slingshot 输出 | ✅ 已完成，保留作基础结果 |
| `fig05_cellchat_prepare.R` | epn2_r | Figure 5A + CellChat 对象（稀疏合并 + CellChat 推断） | ✅ 已完成，输出 `output/fig5/cellchat_merged.rds` |
| `fig05_cellchat_downstream.R` | epn2_r | Figure 5B/C + Supp Fig 8（CellChat 下游全量图 + 可读版图） | ✅ 已整理，输出 `output/fig5_fullsize/` 和 `output/fig5_readable/` |
| `fig06_pyscenic_run.py` / `fig06_pyscenic_downstream.py` | epn2 | Figure 6（regulon heatmap） | ✅ 已完成 GRN、ctx、aucell 和下游热图，输出 `output/fig6/` |
| `fig04_survival_template.R` | epn2_r | Figure 4C（通用 KM 模板） | 🟡 保留外部队列接口 |
| `fig04_gillen_survival_km.R` | epn2_r | Figure 4C（Gillen GSE125861 OS/PFS KM） | ✅ 已完成，输出 OS/PFS trajectory High/Low KM |
| `fig04_gillen_scrna_validation.py` | epn2 | Figure 4B/4D（Gillen GSE125969 external scRNA validation） | 🟡 已完成 cell composition + trajectory score primary/recurrent 验证；recurrent n=3 且全为 ST-RELA |
| `fig04_gillen_scrna_pseudobulk_lr.py` | epn2 | Figure 4D / 5D-E（GSE125969 pseudobulk DEG + key LR expression） | 🟡 已完成 all-subtype pseudobulk DEG 和 LR expression heatmap；strict subtype-adjusted 结论受 recurrent n=3 限制 |
| `fig04_gillen_scrna_deg_go.R` | epn2_r | Figure 4D（GSE125969 recurrent-high GO） | 🟡 已完成 recurrent-high pseudobulk GO BP 富集 |
| `fig04_gojo_scrna_validation.py` | epn2 | Figure 4B/4D/5D-E（Gojo GSE141460 external scRNA validation） | 🟡 已完成 composition、trajectory、pseudobulk DEG 和 LR expression proxy；trajectory recurrent 中位数更高但不显著 |
| `fig04_gojo_scrna_deg_go.R` | epn2_r | Figure 4D（Gojo GSE141460 PF-A recurrent-high GO） | 🟡 PF-A recurrent-high GO 未过富集阈值，保留阴性/保守结果 |
| `fig04_external_scrna_meta_validation.py` | epn2 | Figure 4D/5D-E（Gillen+Gojo low-memory external scRNA meta-validation） | 🟡 已完成 sample-level trajectory model、DEG direction concordance 和 key LR direction concordance；不重新读取原始大矩阵，适合本地 11 GB RAM |
| `fig02_infercnv_build_gene_order.py` / `fig02_infercnv_run_gte009.R` / `fig02_infercnv_genelevel_panels.R` | epn2 / infercnv_r / epn2_r | Figure 2D-F（GTE009 gene-level inferCNV strict CNV 候选） | ✅ WHU HPC inferCNV 已完成并拉回 final object；gene-level panel 脚本可从 final inferCNV object 生成 DYNC2H1 volcano、chr11 group heatmap、chr11 delta profile、sampled-cell heatmap、组合 panel 和统计表 |
| — | — | Figure 2B RNA velocity | 缺 BAM，跳过 |
| — | — | Supp. Fig 2 WES CNV | 缺 WES，跳过 |

## 4. 路径约定

参见仓库根 `.github/copilot-instructions.md` 和 `workspace_paths.py` / `workspace_paths.R`：

- Python：`from workspace_paths import raw_data_path, output_path`
- R：`source("workspace_paths.R")`，使用 `raw_data_path()` / `output_path()`

## 5. 已知坑

- pip 26.0.1 在 Python 3.10 上有 `Link.from_json` bug → 必须降到 `pip<25`
- velocyto 0.17.17 不支持 PEP 517 build isolation，必须 `pip install --no-build-isolation velocyto`，且事先装好 `numpy<2` + `cython`
- velocyto 需要 numpy<2（不兼容 numpy 2.x）
- Monocle2 已不在 CRAN，需从 Bioconductor 安装，且要 R ≤ 4.2 才编译顺利
- 当前 R 4.4 环境中 Monocle2 与新版 dplyr/igraph 存在实际运行兼容问题（`group_by_`/`select_`/`nei` 等旧接口），本项目最终移除 Monocle2 运行路线，Figure 2C 改用 Slingshot 复现核心拟时序结论
- inferCNV 在 BiocManager 上，依赖较重
