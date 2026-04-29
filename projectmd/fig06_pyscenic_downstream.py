#!/usr/bin/env python
"""fig06_pyscenic_downstream.py — Fig 6: regulon AUC 热图 + 细胞型特异 regulon
输入: output/fig6/auc_mtx.loom (来自 fig06_pyscenic_run.py 的 aucell 结果)
输出: output/fig6/Fig6_regulon_AUC_heatmap.pdf
       output/fig6/regulon_celltype_specific.csv
"""
import sys, os
from pathlib import Path
import numpy as np, pandas as pd
import h5py
import matplotlib.pyplot as plt, seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workspace_paths import output_path

FIG6 = Path(output_path('fig6'))
auc_loom = FIG6 / 'auc_mtx.loom'
meta_csv = FIG6 / 'scenic_input' / 'meta.csv'

if not auc_loom.exists():
    print(f'缺 {auc_loom}，请先跑完 fig06_pyscenic_run.py'); sys.exit(1)


def decode_values(values):
    return [v.decode('utf-8') if isinstance(v, (bytes, bytearray)) else str(v) for v in values]


# pySCENIC 的 aucell loom 将 regulon 活性保存在 col_attrs/RegulonsAUC，
# 主矩阵仍是 gene × cell 表达矩阵，不能当成 AUC 读取。
with h5py.File(str(auc_loom), 'r') as h5:
    if 'RegulonsAUC' not in h5['col_attrs']:
        print(f'缺 col_attrs/RegulonsAUC: {auc_loom}'); sys.exit(1)
    regulons_auc = h5['col_attrs']['RegulonsAUC'][:]
    regulon_names = list(regulons_auc.dtype.names or [])
    cell_ids = decode_values(h5['col_attrs']['CellID'][:])
    if not regulon_names:
        print(f'RegulonsAUC 没有 regulon 字段: {auc_loom}'); sys.exit(1)
    auc_cells = pd.DataFrame(
        {name: regulons_auc[name] for name in regulon_names},
        index=cell_ids,
        dtype=np.float32,
    )
    auc = auc_cells.T
print(f'AUC: {auc.shape}')

# meta
meta = pd.read_csv(meta_csv, index_col=0)
common_cells = auc.columns.intersection(meta.index)
auc = auc.loc[:, common_cells]
meta = meta.loc[common_cells]
ct_col = next((c for c in ['Brief_cluster','celltype','cluster'] if c in meta.columns), meta.columns[0])
print(f'用 {ct_col} 分组 ({meta[ct_col].nunique()} 类)')

# 每细胞型平均 AUC
mean_auc = auc.T.groupby(meta[ct_col]).mean().T   # regulon × celltype
mean_auc.to_csv(FIG6 / 'regulon_AUC_mean_by_celltype.csv')

# Top 10 specific regulons per celltype (Z-score)
z = mean_auc.sub(mean_auc.mean(axis=1), axis=0).div(mean_auc.std(axis=1)+1e-6, axis=0)
top = {ct: z.nlargest(10, ct).index.tolist() for ct in z.columns}
pd.DataFrame.from_dict(top, orient='index').T.to_csv(FIG6 / 'regulon_celltype_top10.csv')

# Heatmap of top regulons across celltypes
top_regs = sorted(set(g for v in top.values() for g in v))
plt.figure(figsize=(max(6, mean_auc.shape[1]*0.5), max(8, len(top_regs)*0.18)))
sns.heatmap(z.loc[top_regs], cmap='RdBu_r', center=0, cbar_kws={'label':'Z-score AUC'})
plt.title('Cell-type specific regulons (SCENIC AUC)')
plt.tight_layout()
plt.savefig(FIG6 / 'Fig6_regulon_AUC_heatmap.pdf')
plt.close()

print(f'✓ Fig 6 完成 → {FIG6}/Fig6_regulon_AUC_heatmap.pdf')
