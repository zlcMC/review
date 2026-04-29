# Figure 6 pySCENIC regulon 分析：从矩阵开始理解

这一步对应论文中 SCENIC / regulon 分析。原文说法可以拆成一句话：把 4 个 EPN 样本整合后的 log-normalized 表达矩阵输入 pySCENIC，推断每种细胞里哪些转录因子调控程序更活跃。

在本项目里，相关文件是：

| 目的 | 文件 |
|---|---|
| 主流程脚本 | `projectmd/fig06_pyscenic_run.py` |
| 下游热图脚本 | `projectmd/fig06_pyscenic_downstream.py` |
| 输入矩阵目录 | `output/fig6/scenic_input/` |
| pySCENIC 输出目录 | `output/fig6/` |
| 最终热图 | `output/fig6/Fig6_regulon_AUC_heatmap.pdf` |

## 1. 先把这一步放回单细胞分析主线

前面几步大致是：

| 分析层次 | 你在问的问题 | 矩阵是什么 |
|---|---|---|
| Figure 1 聚类/注释 | 这个细胞像哪类细胞？ | 基因 x 细胞表达矩阵 |
| Figure 2/3 拟时序 | 肿瘤细胞从哪种状态走向哪种状态？ | 细胞 x 轨迹分数 / pseudotime |
| Figure 5 CellChat | 哪些细胞通过配体-受体互相说话？ | 发送细胞 x 接收细胞 x 通路 |
| Figure 6 pySCENIC | 哪些转录因子调控程序在某类细胞中活跃？ | regulon x 细胞 / 细胞类型 |

所以 Figure 6 这一步不是重新聚类，也不是找 marker gene，而是从表达矩阵中进一步问：这些基因表达变化背后，可能是哪批转录因子程序在驱动？

## 2. 一个最小例子：表达矩阵长什么样

假设只有 5 个基因、4 个细胞：

| gene / cell | Cell_1 NSC | Cell_2 NSC | Cell_3 Mic | Cell_4 EC |
|---|---:|---:|---:|---:|
| SOX2 | 5.1 | 4.8 | 0.1 | 0.0 |
| HES5 | 3.2 | 3.5 | 0.0 | 0.1 |
| C1QA | 0.0 | 0.1 | 5.8 | 0.0 |
| PECAM1 | 0.0 | 0.0 | 0.1 | 4.9 |
| JUN | 1.8 | 2.1 | 2.4 | 1.7 |

单细胞数据最底层就是这种矩阵。行是基因，列是细胞，格子里的数字是这个基因在这个细胞里的表达量。真实项目里大很多：当前 `output/fig6/scenic_input/genes.txt` 有 16,448 个基因，`cells.txt` 有 35,102 个细胞。

也就是：

| 真实输入矩阵 | 数量 |
|---|---:|
| genes | 16,448 |
| cells | 35,102 |
| 理论格子数 | 16,448 x 35,102 |

这个矩阵被保存为 `output/fig6/scenic_input/integrated.loom`。`loom` 可以理解成专门给单细胞大矩阵用的文件盒子：主矩阵 + 行注释 + 列注释。

## 3. pySCENIC 三步分别在做什么

`projectmd/fig06_pyscenic_run.py` 跑的是标准三步：

```text
表达矩阵
  -> grn
  -> ctx
  -> aucell
  -> regulon 活性矩阵
```

### 3.1 grn：从表达相关性推 TF-target 候选边

输入：表达矩阵 + 转录因子列表。

输出：`output/fig6/adj.tsv`。

这个文件是边表，不是传统二维矩阵。它长这样：

| TF | target | importance |
|---|---|---:|
| JUNB | ZFP36 | 59.04 |
| ASCL1 | HES6 | 45.53 |
| EGR1 | FOS | 45.52 |
| FOS | JUNB | 38.76 |

可以把它想成一个稀疏矩阵：

| TF / target | ZFP36 | HES6 | FOS | JUNB |
|---|---:|---:|---:|---:|
| JUNB | 59.04 | 0 | 0 | 0 |
| ASCL1 | 0 | 45.53 | 0 | 0 |
| EGR1 | 0 | 0 | 45.52 | 0 |
| FOS | 0 | 0 | 0 | 38.76 |

`importance` 不是 p 值，可以先理解为模型认为 TF 和 target 之间联系强不强。它只说明表达模式上像有关联，还不能证明这个 TF 真的直接调控了这个基因。

当前真实结果中，`adj.tsv` 大约有 1,487,125 行，说明第一步先产生了很多候选 TF-target 关系。

### 3.2 ctx：用 motif 证据筛选候选边，形成 regulon

输入：`adj.tsv` + cisTarget motif ranking 数据库 + motif 注释。

输出：`output/fig6/regulons.csv`。

这一步回答的是：某个 TF 的候选 target 附近，有没有这个 TF 对应的 DNA motif 富集？如果有，调控关系更可信。

regulon 可以这样理解：

```text
一个 regulon = 一个转录因子 + 一组被它调控的 target genes
```

例子：

| regulon | target genes 的概念 |
|---|---|
| BRCA1(+) | KPNB1, GGCT, CDC6, MCM5, PCNA, ... |
| SOX2(+) | 一批可能受 SOX2 调控的基因 |
| MEF2C(+) | 一批可能受 MEF2C 调控的基因 |

名字里的 `(+)` 表示 activating regulon，也就是这个转录因子活跃时，target genes 整体倾向于更高表达。

### 3.3 aucell：给每个细胞打 regulon 活性分

输入：原始表达矩阵 + regulon 列表。

输出：`output/fig6/auc_mtx.loom`。

这一步生成最关键的矩阵：regulon x cell。

一个小例子：

| regulon / cell | Cell_1 NSC | Cell_2 NSC | Cell_3 Mic | Cell_4 EC |
|---|---:|---:|---:|---:|
| SOX2(+) | 0.82 | 0.79 | 0.05 | 0.02 |
| MEF2C(+) | 0.10 | 0.12 | 0.88 | 0.20 |
| BCL6B(+) | 0.03 | 0.04 | 0.08 | 0.91 |

这张表的含义不是“SOX2 基因表达量”，而是“SOX2 regulon 这套 target genes 在这个细胞里整体是否排在表达谱前列”。

临床直觉类比：单个 marker 像一次化验指标，regulon AUC 更像一组指标组成的综合评分。

## 4. 为什么下游脚本要再按细胞类型求平均

`auc_mtx.loom` 是 regulon x 35,102 个细胞。论文图通常不会直接画 35,102 列，因为读不动，所以 `projectmd/fig06_pyscenic_downstream.py` 做了两件事：

1. 把每个细胞的 regulon AUC 读出来。
2. 按 metadata 里的 `Brief_cluster` 分组，求每类细胞的平均 AUC。

于是得到 `output/fig6/regulon_AUC_mean_by_celltype.csv`：regulon x cell type。

真实输出的开头类似：

| regulon | Ast | EC | Epe | Mal | Mic | NSC | Neu | OPC | Oli | Per | TC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ATF3(+) | 0.072 | 0.061 | 0.042 | 0.057 | 0.104 | 0.076 | 0.045 | 0.051 | 0.039 | 0.062 | 0.075 |
| BCL6B(+) | 0.0046 | 0.202 | 0.0037 | 0.0057 | 0.010 | 0.011 | 0.0066 | 0.012 | 0.0040 | 0.057 | 0.035 |
| CEBPA(+) | 0.034 | 0.038 | 0.024 | 0.019 | 0.216 | 0.017 | 0.022 | 0.022 | 0.035 | 0.064 | 0.082 |

怎么读？

- `BCL6B(+)` 在 EC 列是 0.202，明显高于很多其他细胞类型，说明这个 regulon 在 EC 中更活跃。
- `CEBPA(+)` 在 Mic 列是 0.216，提示它可能是 Mic 相关调控程序。
- 这些是活性分，不是差异表达 p 值；后续解释要结合细胞身份、target genes 和文献背景。

## 5. top10 文件怎么读

`output/fig6/regulon_celltype_top10.csv` 是每个细胞类型挑出最特异的 10 个 regulon。

真实输出开头：

| rank | Ast | EC | Epe | Mal | Mic | NSC |
|---:|---|---|---|---|---|---|
| 0 | MAFG(+) | BCL6B(+) | RELA(+) | POU2F1(+) | MEF2C(+) | MEIS1(+) |
| 1 | E2F8(+) | SOX18(+) | STAT5A(+) | SOX3(+) | MAF(+) | EGR1(+) |
| 2 | JUN(+) | FLI1(+) | SMC3(+) | STAT5A(+) | IRF5(+) | HES5(+) |

这张表适合当“索引”：你想看 NSC-like 细胞，就先看 NSC 列；想看 malignant cells，就看 Mal 列。

## 6. 和论文问题怎么连起来

论文不只是想画一张 regulon 热图，而是想把调控网络和复发相关的细胞通讯结合起来。原文大意是：

1. 在 4 个整合 EPN 样本上跑 SCENIC。
2. 找出 NSC-like cells 中活性显著高的 regulons。
3. 再结合 CellChat 中 NSC-like 相关的 ligand/receptor。
4. 看这些基因是否共同落在相似的 GO/KEGG 生物学通路里。

所以你可以把 Figure 6 pySCENIC 理解成 Figure 5 CellChat 之后的“调控层补充”：

| Figure 5 | Figure 6 / Supplementary regulon |
|---|---|
| 细胞之间怎么交流 | 细胞内部哪些 TF 程序活跃 |
| ligand-receptor | TF-target regulon |
| sender x receiver | regulon x cell type |
| 偏细胞外信号 | 偏细胞内转录调控 |

## 7. 当前项目的运行状态

当前 `output/fig6/` 已经包含主结果：

| 文件 | 当前含义 |
|---|---|
| `adj.tsv` | grn 得到的 TF-target 候选边表，约 148 万行 |
| `regulons.csv` | ctx 得到的 motif 支持 regulon，约 6,298 行 motif/regulon 记录 |
| `auc_mtx.loom` | aucell 得到的 regulon 活性结果 |
| `regulon_AUC_mean_by_celltype.csv` | regulon x cell type 的均值矩阵，约 92 个 regulon |
| `regulon_celltype_top10.csv` | 每类细胞 top regulon |
| `Fig6_regulon_AUC_heatmap.pdf` | 细胞类型特异 regulon 热图 |

如果重新运行，命令是：

```bash
conda run -n epn2 python projectmd/fig06_pyscenic_run.py
conda run -n epn2 python projectmd/fig06_pyscenic_downstream.py
```

## 8. 初学者最容易混淆的三个点

### 8.1 regulon 不是单个基因

`SOX2` 是一个基因 / 转录因子，`SOX2(+)` 是一个 regulon，代表 SOX2 加上一组 target genes。看 regulon AUC 时，不是在看 SOX2 自己表达多高，而是在看 SOX2 这套调控程序整体有多活跃。

### 8.2 AUC 不是表达量

AUC 是一种排序型分数。它问的是：某个 regulon 的 target genes 是否集中出现在这个细胞表达谱的高排名区域。这个思路比单个基因更稳，因为它利用了一组 target genes。

### 8.3 热图颜色通常来自标准化后的矩阵

`fig06_pyscenic_downstream.py` 里先得到平均 AUC，再对每个 regulon 跨细胞类型做 Z-score：

```text
z = (某 regulon 在某细胞类型的平均 AUC - 该 regulon 在所有细胞类型的均值) / 标准差
```

所以热图颜色表示“这个 regulon 相对其他细胞类型是不是更高”，不是原始 AUC 的绝对值。

## 9. 推荐你学习时的阅读顺序

1. 先看 `output/fig6/regulon_AUC_mean_by_celltype.csv`，理解 regulon x cell type 矩阵。
2. 再看 `output/fig6/regulon_celltype_top10.csv`，找每类细胞的代表 regulon。
3. 再回头看 `output/fig6/adj.tsv`，理解 TF-target 边表。
4. 最后看 `output/fig6/regulons.csv`，理解 motif 证据如何把候选边变成 regulon。
5. 打开 `output/fig6/Fig6_regulon_AUC_heatmap.pdf`，把热图颜色和前面的矩阵对应起来。
