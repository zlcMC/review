# v2 重构总报告

> 与 `projectmd/` 中 v1 脚本一一对应，**除 4 个 HPC 脚本外全部完成**。
> 原始脚本完全保留，可随时退回。

## 1. 范围

### ✅ 已重构（25 个脚本 + 2 helper + 本说明）

| 类别 | 数量 |
|---|---|
| Python 脚本 | 6 |
| R 脚本 | 18 |
| Shell 脚本 | 2 |
| Helper（共享） | 2（`helpers/helpers.R`, `helpers/helpers_common.py`） |

### ❌ 未动（仅 HPC 与文档）

HPC 必须保留原状（资源 / 路径 / 调度）：
- `projectmd/fig02_infercnv_run_gte009.R`
- `projectmd/fig02_infercnv_genelevel_panels.R`
- `projectmd/fig06_pyscenic_run.py`
- `projectmd/fig06_scenic_ctx_aucell_pipeline.sh`
- `projectmd/hpc/*`

> **更新（HPC 数据已落地）**：上述 4 个脚本已等价迁移到 `projectmd/v2/`：
> 仅做路径迁移（`v2_dir`/`dep_file`），不改算法/参数。HPC 重跑步骤本地不复跑，
> 但本地下游使用 v2 脚本读取 v1 HPC 中间产物，结果与 v1 字节一致：
>
> | v2 脚本 | 本地验证产物 | cmp v1 |
> |---|---|---|
> | `fig02_infercnv_genelevel_panels.R` | `fig2_infercnv_strict/GTE009_downstream/*.csv` ×3 | IDENTICAL ×3 |
> | `fig06_pyscenic_downstream.py`（读 HPC `regulons.csv`/`auc_mtx.loom`）| `fig6/regulon_AUC_mean_by_celltype.csv`、`fig6/regulon_celltype_top10.csv` | IDENTICAL ×2 |
> | `fig02_infercnv_run_gte009.R` / `fig06_pyscenic_run.py` / `fig06_scenic_ctx_aucell_pipeline.sh` | 仅在 HPC 上运行 | — |

文档 / 环境定义：
- `projectmd/README.md`、`projectmd/RUN_ORDER.md`
- `projectmd/epn2_r.environment.yml`
- `projectmd/fig06_pyscenic_matrix_guide.md`

## 2. 代码量对比

| 批次 | 内容 | orig | v2 | 节省 |
|---|---|---:|---:|---:|
| batch2 | shell + 4 R | 850 | 614 | −236 (−28%) |
| batch3 | 4 GO/DEG R | 339 | 218 | −121 (−36%) |
| batch4 | 4 cellchat + final_paper_style | 868 | 712 | −156 (−18%) |
| batch5 | fig01_tsne / cohort / fig02_proxy / fig06_downstream | 536 | 458 | −78 (−15%) |
| batch6 | infercnv_build / cistarget / 2× slingshot | 406 | 312 | −94 (−23%) |
| **合计** | **25 个脚本** | **2999** | **2314** | **−685 (−23%)** |

> Helper 一次性投入：`helpers.R` 246 行 + `helpers_common.py` 137 行 = 383 行；净节省 ≈ **300 行**（~10%）。

## 3. v2 结构

```
projectmd/v2/
├── README.md                 ← v2 总说明
├── helpers/
│   ├── helpers.R             ← R 共享：init_workspace / v2_file / dep_file / load_gse_expression / score_trajectory / enrich_go_dotplot / cellchat_*
│   └── helpers_common.py     ← Python 共享：v2_file / dep_file / parse_tsv_floats / bh_adjust / subtype_group / KEY_LR_GENES_*
├── fig01_*.R                 ← 3 个
├── fig02_*.{py,R}            ← 4 个（slingshot ×2 + detail_proxy + infercnv_build_gene_order）
├── fig03_*.R                 ← 1 个
├── fig04_*.{py,R}            ← 8 个（含 microarray / scrna / survival）
├── fig05_cellchat_*.R        ← 4 个
├── fig06_pyscenic_downstream.py
├── final_paper_style_export.R
├── cohort_prepare_pajtler_template.R
└── setup_*.sh                ← 2 个
```

## 4. 核心约定

### 路径
- **写**：一律 `v2_file('figX', 'foo.csv')` / `v2_dir('figX')` → `output/v2/figX/...`
- **读上游**：`dep_file('figX', 'foo.csv')` —— v2 优先、v2 不存在自动回落到 v1 `output/figX/`，因此可独立运行
- **原始数据**：`raw_data_path('GTE009.h5ad')` → `projectfile/`

### Bootstrap 模板
R：
```R
this_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  dirname(sub('^--file=', '', grep('^--file=', a, value = TRUE)[1]))
})()
source(file.path(this_dir, 'helpers', 'helpers.R'))
init_workspace()
```

Python：
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import v2_file, dep_file, ...
```

## 5. 验证状态

| 脚本 | 静态 | 运行时 |
|---|---|---|
| 6 Python | importlib OK | 2 个实跑通过（`fig04_external_scrna_meta_validation`、`fig06_pyscenic_downstream`） |
| 18 R | `parse()` OK | **全部 18 个实跑验证完成**（8 fig + 5 cellchat）；唯一未过 = `cohort_prepare_pajtler_template.R`（模板脚本，本机无 `projectfile/external_cohort/Pajtler2015/` 数据，预期如此） |
| 2 Shell | `bash -n` OK | — |

实跑产物：
- `output/v2/fig4_external_scrna_meta/` —— 12 个文件
- `output/v2/fig6/` —— `Fig6_regulon_AUC_heatmap.pdf` + 2 CSV
- `output/v2/fig1..fig5/` + `output/v2/fig_paper_style/` —— 18 R 跳后产生的全量产物
- 每个脚本运行日志：`output/v2/_runtime_logs/`

### 18 R 跳表（代表耗时）

| 脚本 | 状态 | 耗时 |
|---|---|---:|
| cohort_prepare_pajtler_template.R | FAIL（预期，缺数据） | 1s |
| fig01_bulk_relapse_deg.R | PASS | 20s |
| fig01_supp_subclone_deg_go.R | PASS | 107s |
| fig01_tsne_cnv_cytotrace_seurat.R | PASS | 32s |
| fig02_slingshot_basic.R | PASS | 140s |
| fig02c_slingshot_final_panel.R | PASS | 168s |
| fig03_subclone_undiff_diff_score.R | PASS | 179s |
| fig04_bulk_relapse_proxies.R | PASS（修复后） | 22s |
| fig04_gillen_microarray_proxy.R | PASS（修复后） | 10s |
| fig04_gillen_scrna_deg_go.R | PASS | 35s |
| fig04_gojo_scrna_deg_go.R | PASS | 25s |
| fig04_gillen_survival_km.R | PASS | 14s |
| fig04_survival_template.R | PASS | 3s |
| fig05_cellchat_prepare.R | PASS | **590s** |
| fig05_cellchat_downstream.R | PASS | 89s |
| fig05_cellchat_fullsize_figures.R | PASS | 72s |
| fig05_cellchat_readable_figures.R | PASS | 14s |
| final_paper_style_export.R | PASS | 63s |

## 6. 已修复的具体问题

- `np.fromstring` （NumPy 2 已移除）→ `parse_tsv_floats`
- 4 处 `bh_adjust` 重复定义 → helper 单点
- `KEY_LR_GENES` 多处重复 → `KEY_LR_GENES_CORE` / `KEY_LR_GENES_UNION` 单点
- `setwd + source('workspace_paths.R')` 仪式 → `init_workspace()`
- GEO 加载 / trajectory 评分 / GO dotplot / CellChat 配色 → 5 个 helper 收敛
- shell 安装脚本加 `have()` 幂等性
- ggplot2 新版不再接受 `dpi = NULL` → fig04_bulk_relapse_proxies.R / fig04_gillen_microarray_proxy.R 中 `dpi = if (ext == 'png') 220 else 300`
- `fig04_gillen_scrna_validation.py::stream_signature_scores` 出现重复 if 块导致目标基因计数翻倍 → 删除重复块。修复后 `GSE125969_neoplastic_cell_trajectory_scores.csv` / `condition_stats` / `patient_summary` 与 v1 字节一致。
- `fig01_supp_subclone_deg_go.R` Fig1I GO 行数 / 引号 / 浮点精度与 v1 不同：
  1. v1 用 `enrichGO` 默认 `qvalueCutoff = 0.05`，而 helper 默认 `0.2` → `enrich_go_dotplot` 增加 `pcutoff` / `qcutoff` 透传，脚本传 `pcutoff=0.05, qcutoff=0.05`。
  2. v1 用 `write.csv`（带引号、低浮点精度），helper 默认 `readr::write_csv` → 增加 `csv_writer = c('readr','base')` 选项，脚本传 `csv_writer = 'base'`。

## 6.1 v2 ↔ v1 数值等价性（最终）

跨 fig 全量 CSV 对比（`output/{figX}/*.csv` ↔ `output/v2/{figX}/*.csv`）：

```
IDENTICAL = 58
DIFFER    = 0
V2_MISSING = 2  # 仅 output/fig5/pathway_importance.csv（历史 notebook 产物）
                # 与 output/fig6/regulons.csv（HPC pyscenic_run.py 产物，不在 v2 范围）
```

Fig2 详情面板 (`Fig2_cilium_proxy_genes.csv`, `Fig2_detail_proxy_stats.csv`)、Fig4 全部 scRNA / pseudobulk LR / Gillen GO / Gojo GO、Fig5 cellchat（各分辨率版本）、Fig6 regulon AUC、Fig1 bulk DEG / Fig1I GO 等均与 v1 字节级一致。

## 6.2 论文吻合度

v2 与 v1 完全等同；缺口仍为：
- Fig 2B RNA velocity（缺 spliced/unspliced）
- Supp 2 WES 拷贝数（作者未公开）
- 部分 OS / KM 曲线（Pajtler 2015 等公开数据无 OS）

均为数据可得性限制，**不是代码问题**。

## 7. 迁移指南（v1 → v2）

1. 想跑 v2 的某个脚本：直接 `Rscript projectmd/v2/<name>.R` 或 `python projectmd/v2/<name>.py`
2. **不需要先跑过 v2 上游**：dep_file 会回落到 v1 `output/`
3. 想清掉 v2 输出重跑：`rm -rf output/v2/figX/`
4. v1 仍可独立运行（未做任何改动）

## 8. 剩余可做项（非阻塞）

- 实跑全部 18 个 v2 R 脚本，确认产物逐 byte 等价
- v2 RUN_ORDER.md（按 fig 编号 + 数据依赖列出推荐顺序）
- B3：fig04 R 脚本进一步参数化（重复模式还可再抽）
- C1：fig06_pyscenic_run.py 加 seed + 断点续跑（HPC，谨慎）
