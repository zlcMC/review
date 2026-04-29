#!/usr/bin/env python
"""fig06_pyscenic_run.py (v2) — pySCENIC GRN→ctx→AUCell.
v1 等效；切到 v2_dir/dep_file 路径与 helpers_common；输出 output/v2/fig6/。
HPC 级耗时（GRN 数小时），本地不重跑。
"""
import os, sys, subprocess, textwrap
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / 'helpers'))
from helpers_common import v2_dir, dep_file  # noqa: E402

FIG6      = v2_dir('fig6')
CIS       = v2_dir('fig6', 'cistarget')
SCENIC_IN = v2_dir('fig6', 'scenic_input')

NEEDED = {
    'rankings': 'hg19-500bp-upstream-7species.mc9nr.feather',
    'motifs':   'motifs-v9-nr.hgnc-m0.001-o0.0.tbl',
    'tfs':      'allTFs_hg38.txt',
}

def resolve_db(name: str) -> Path:
    """v2 优先，否则回退 v1 cistarget 目录。"""
    p_v2 = CIS / name
    if p_v2.exists():
        return p_v2
    p_v1 = Path(dep_file('fig6', 'cistarget', name))
    return p_v1 if p_v1.exists() else p_v2

dbs = {k: resolve_db(v) for k, v in NEEDED.items()}
missing = [str(p) for p in dbs.values() if not p.exists()]
if missing:
    print(textwrap.dedent(f"""
    ===========================================================
    缺少 cisTarget 数据库文件（放到 {CIS}/ 或 v1 路径）：
      {chr(10).join('  - ' + m for m in missing)}
    下载命令：
      wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather
      wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
      wget https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/allTFs_hg38.txt
    ===========================================================
    """).strip()); sys.exit(1)

# ---- 1) 载入 fig05_cellchat_prepare 输出 → loom ----
loom_file = SCENIC_IN / 'integrated.loom'
if not loom_file.exists():
    mtx_rds = Path(dep_file('fig6', 'scenic_input', 'integrated_logmat.rds'))
    if not mtx_rds.exists():
        # 也允许 v1 / v2 cellchat prepare 路径
        mtx_rds = SCENIC_IN / 'integrated_logmat.rds'
    if not mtx_rds.exists():
        print(f'缺少 {mtx_rds}，请先跑 fig05_cellchat_prepare.R'); sys.exit(1)
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

    import scipy.io, numpy as np, loompy
    mtx = scipy.io.mmread(str(SCENIC_IN / 'expr.mtx')).tocsc()
    genes = [l.strip() for l in open(SCENIC_IN / 'genes.txt')]
    cells = [l.strip() for l in open(SCENIC_IN / 'cells.txt')]
    print(f'矩阵: {mtx.shape}; 分块写 loom ...')
    chunk = 2000
    first = mtx[:, :chunk].toarray().astype(np.float32)
    loompy.create(str(loom_file), first,
                  {'Gene': np.array(genes)},
                  {'CellID': np.array(cells[:chunk])})
    del first
    with loompy.connect(str(loom_file)) as ds:
        for i in range(chunk, mtx.shape[1], chunk):
            j = min(i + chunk, mtx.shape[1])
            block = mtx[:, i:j].toarray().astype(np.float32)
            ds.add_columns(block, col_attrs={'CellID': np.array(cells[i:j])})
            del block
            print(f'  写入 {j}/{mtx.shape[1]}')
    print('loom:', loom_file)

# ---- 2) GRN / ctx / AUCell ----
adj_tsv  = FIG6 / 'adj.tsv'
reg_csv  = FIG6 / 'regulons.csv'
auc_loom = FIG6 / 'auc_mtx.loom'

def sh(cmd):
    print('$', ' '.join(cmd)); subprocess.check_call(cmd)

if not adj_tsv.exists():
    sh(['pyscenic', 'grn', str(loom_file), str(dbs['tfs']),
        '-o', str(adj_tsv), '--num_workers', '4'])
if not reg_csv.exists():
    sh(['pyscenic', 'ctx', str(adj_tsv), str(dbs['rankings']),
        '--annotations_fname', str(dbs['motifs']),
        '--expression_mtx_fname', str(loom_file),
        '-o', str(reg_csv), '--num_workers', '4'])
if not auc_loom.exists():
    sh(['pyscenic', 'aucell', str(loom_file), str(reg_csv),
        '-o', str(auc_loom), '--num_workers', '4'])

print('\n✓ pySCENIC 完成:', auc_loom)
