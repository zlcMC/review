#!/usr/bin/env python
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python (epn2)
#     language: python
#     name: epn2
# ---

# %%
"""fig06_pyscenic_run.py — 对应论文 Figure 6（pySCENIC regulon 分析）

论文 Methods：
  - 在 4 个整合样本上运行 pySCENIC 默认流程
  - 输入：log-normalized 表达矩阵（由 fig05_cellchat_prepare.R 导出）
  - grn → ctx (cisTarget motif ranking) → aucell

运行：conda run -n epn2 python projectmd/fig06_pyscenic_run.py

注意：
  1) 需要先下载 cisTarget 数据库 (hg19/hg38) 和 motif annotation (约 5GB)。
     脚本会检查 output/fig6/cistarget/ 是否存在，不存在就打印下载指令而不是自动下。
  2) 需要 pyarrow / pyscenic。已在 epn2 环境中 (pySCENIC 0.12.1)。
"""
import os, sys, subprocess, textwrap
from pathlib import Path
sys.path.insert(0, os.getcwd())
from workspace_paths import output_path

# %%
FIG6 = Path(output_path('fig6'))
CIS = FIG6 / 'cistarget'
SCENIC_IN = FIG6 / 'scenic_input'
CIS.mkdir(parents=True, exist_ok=True)
SCENIC_IN.mkdir(parents=True, exist_ok=True)

# %%
# ---- 0) 检查 cisTarget 数据库 ----
# 论文用 hg19 (CellRanger hg19-3.0.0)，选用 hg19 数据库最匹配
needed = {
    'rankings': 'hg19-500bp-upstream-7species.mc9nr.feather',
    'motifs':   'motifs-v9-nr.hgnc-m0.001-o0.0.tbl',
    'tfs':      'allTFs_hg38.txt',  # pySCENIC 自带 TF list，hg38/hg19 同名
}
missing = [v for v in needed.values() if not (CIS / v).exists()]
if missing:
    msg = textwrap.dedent(f"""
    ===========================================================
    缺少 cisTarget 数据库文件 (放到 {CIS}/)：
      {chr(10).join('  - ' + m for m in missing)}

    下载命令（约 1GB ranking + 5MB motifs + 2KB TFs）：
      cd {CIS}
      wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather
      wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
      wget https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/allTFs_hg38.txt

    下载完成后再运行本脚本。
    ===========================================================
    """).strip()
    print(msg)
    sys.exit(1)

# %%
# ---- 1) 载入 Fig 5 prepare 步骤导出的矩阵 → loom ----
loom_file = SCENIC_IN / 'integrated.loom'
if not loom_file.exists():
    print('转换 R 矩阵 → loom ...')
    # 通过 anndata2ri / rpy2 太重；改用 R 端导出 csv+meta，再 loompy 组装
    mtx_rds = SCENIC_IN / 'integrated_logmat.rds'
    if not mtx_rds.exists():
        print(f'缺少 {mtx_rds}，请先跑 fig05_cellchat_prepare.R'); sys.exit(1)
    # 用 R 一行导出为 10X mtx 再读
    r_export = SCENIC_IN / '_export_mtx.R'
    r_export.write_text(textwrap.dedent(f"""
        setwd('{os.getcwd()}')
        source('workspace_paths.R')
        x <- readRDS('{mtx_rds}')
        Matrix::writeMM(x$expr, '{SCENIC_IN}/expr.mtx')
        writeLines(rownames(x$expr), '{SCENIC_IN}/genes.txt')
        writeLines(colnames(x$expr), '{SCENIC_IN}/cells.txt')
        write.csv(x$meta,             '{SCENIC_IN}/meta.csv')
    """).strip())
    subprocess.check_call(['conda', 'run', '-n', 'epn2_r', 'Rscript', str(r_export)])

    import scipy.io, numpy as np, loompy, pandas as pd
    mtx = scipy.io.mmread(str(SCENIC_IN / 'expr.mtx')).tocsc()  # gene × cell sparse
    genes = [l.strip() for l in open(SCENIC_IN / 'genes.txt')]
    cells = [l.strip() for l in open(SCENIC_IN / 'cells.txt')]
    print(f'矩阵: {mtx.shape}; 内存友好分块写 loom ...')
    # pySCENIC 期望 row=Gene, col=Cell —— mtx 已经是 gene × cell
    row_attrs = {'Gene': np.array(genes)}
    col_attrs = {'CellID': np.array(cells)}
    # 分块写避免一次性 dense 化（35K×16K dense = 4.5GB）
    chunk = 2000
    # 先创建一个空 loom（用第一块）
    first = mtx[:, :chunk].toarray().astype(np.float32)
    loompy.create(str(loom_file), first,
                  row_attrs, {'CellID': np.array(cells[:chunk])})
    del first
    with loompy.connect(str(loom_file)) as ds:
        for i in range(chunk, mtx.shape[1], chunk):
            j = min(i + chunk, mtx.shape[1])
            block = mtx[:, i:j].toarray().astype(np.float32)
            ds.add_columns(block, col_attrs={'CellID': np.array(cells[i:j])})
            del block
            print(f'  写入 {j}/{mtx.shape[1]}')
    print(f'loom: {loom_file}')

# %%
# ---- 2) pySCENIC 三步 ----
adj_tsv    = FIG6 / 'adj.tsv'
reg_csv    = FIG6 / 'regulons.csv'
auc_loom   = FIG6 / 'auc_mtx.loom'

# %%
def sh(cmd):
    print('$', ' '.join(cmd)); subprocess.check_call(cmd)

# %%
if not adj_tsv.exists():
    sh(['pyscenic', 'grn', str(loom_file), str(CIS/needed['tfs']),
        '-o', str(adj_tsv), '--num_workers', '4'])
if not reg_csv.exists():
    sh(['pyscenic', 'ctx', str(adj_tsv), str(CIS/needed['rankings']),
        '--annotations_fname', str(CIS/needed['motifs']),
        '--expression_mtx_fname', str(loom_file),
        '-o', str(reg_csv), '--num_workers', '4'])
if not auc_loom.exists():
    sh(['pyscenic', 'aucell', str(loom_file), str(reg_csv),
        '-o', str(auc_loom), '--num_workers', '4'])

# %%
print('\n✓ pySCENIC 完成:', auc_loom)
print('后续作图（regulon activity heatmap 按 Brief_cluster）可在 jupyter 里按需做。')
