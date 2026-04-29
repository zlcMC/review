# 武大超算 (WHU HPC) 运行说明

把吃内存的 pySCENIC (fig6) 步骤搬到武大超算运行，其余步骤继续在本地 WSL 跑。

## 集群信息速览

- 用户名 / 家目录: `zhanglichuan` / `/home/zhanglichuan` (**只有 1 GB 配额**, 别放数据)
- 大盘: `/project/zhanglichuan/` (**500 GB**)
- 登录节点: `swarm01.whu.edu.cn` / `swarm02.whu.edu.cn`
- 文件传输节点: `swarm-xfe.whu.edu.cn` (rsync/scp 走这个更快)
- 校外: 先开武大 SSLVPN
- 调度: SLURM, 我能用的分区 = `hpxg` (16 核 / 90 GB / 24 h 上限)
- `sinfo` 被全局 alias 了, 用 `\sinfo` 跳过

## 工作目录约定 (集群上)

```
/project/zhanglichuan/
├── miniforge3/                      # conda 本体 (装一次)
├── conda_envs/scenic/               # pyscenic 环境
└── wslproject/
    ├── fig6/                        # 输入数据 (从本地 rsync 上来)
    │   ├── scenic_input/            # integrated.loom, expr.mtx, cells.txt, genes.txt
    │   ├── cistarget/               # 3 个数据库文件
    │   └── adj.tsv                  # GRN 输出 (跑完后会产生)
    ├── logs/                        # SLURM 输出
    └── projectmd/hpc/               # SLURM 脚本
```

## 0. 一次性: 同步代码 + 数据到集群

在 **WSL** 里跑一次:

```bash
# 同步本仓库代码 (轻量, 不含 output)
rsync -avzP --exclude='output/' --exclude='projectfile/' --exclude='.venv/' \
    --exclude='.git/' --exclude='read/' \
    /home/zlcmc/wslproject/ \
    zhanglichuan@swarm-xfe.whu.edu.cn:/project/zhanglichuan/wslproject/

# 同步 fig6 输入数据 (本次 ~4 GB)
rsync -avzP \
    /home/zlcmc/wslproject/output/fig6/scenic_input \
    /home/zlcmc/wslproject/output/fig6/cistarget \
    /home/zlcmc/wslproject/output/fig6/adj.tsv \
    zhanglichuan@swarm-xfe.whu.edu.cn:/project/zhanglichuan/wslproject/fig6/
```

## 1. 一次性: 准备 conda 环境

登录 swarm01, 跑:

```bash
cd /project/$USER/wslproject/projectmd/hpc
bash 00_prepare_conda_env.sh
```

会装 miniforge 到 `/project/$USER/miniforge3/`, 并建一个叫 `scenic` 的
conda 环境 (含 pyscenic 0.12.1 + dask<2023). 一次性 装好后永久复用。

## 2. 提交 pySCENIC 作业

```bash
cd /project/$USER/wslproject/projectmd/hpc
sbatch 10_pyscenic_grn.slurm        # 第一步: GRN (最吃内存, ~6-12 h)
# 等它跑完, squeue -u $USER 看进度
sbatch 20_pyscenic_ctx_aucell.slurm # 第二/三步: ctx + aucell (~1-2 h)
```

或者一条龙 (中间不退出):

```bash
JOB1=$(sbatch --parsable 10_pyscenic_grn.slurm)
sbatch --dependency=afterok:$JOB1 20_pyscenic_ctx_aucell.slurm
```

## 3. 把结果拉回本地

跑完后, **WSL** 里:

```bash
rsync -avzP \
    zhanglichuan@swarm-xfe.whu.edu.cn:/project/zhanglichuan/wslproject/fig6/regulons.csv \
    zhanglichuan@swarm-xfe.whu.edu.cn:/project/zhanglichuan/wslproject/fig6/auc_mtx.loom \
    zhanglichuan@swarm-xfe.whu.edu.cn:/project/zhanglichuan/wslproject/fig6/adj.tsv \
    /home/zlcmc/wslproject/output/fig6/
```

之后下游脚本 `projectmd/fig06_pyscenic_downstream.py` 在本地 WSL 跑就行。

## 4. 可选: 提交 inferCNV strict CNV 作业

当前 Fig 6 已经完成；如果现在要把超算挂起来，建议优先挂 Fig 2D-F 的
GTE009 inferCNV strict CNV 候选作业。该任务读取 `projectfile/GTE009_seurat.rds`
里的整数 `RNA/counts` layer，用 `Mic/EC/TC` 作为 reference，`Subclone_1/2`
作为观察组。

先把这个作业额外需要的 RDS 传到集群。注意上面的通用同步命令排除了
`projectfile/`，所以这里必须单独传：

```bash
ssh zhanglichuan@swarm-xfe.whu.edu.cn 'mkdir -p /project/zhanglichuan/wslproject/projectfile'
rsync -avzP \
    /home/zlcmc/wslproject/projectfile/GTE009_seurat.rds \
    zhanglichuan@swarm-xfe.whu.edu.cn:/project/zhanglichuan/wslproject/projectfile/
```

前提：集群上需要有 R + Seurat + inferCNV 环境。推荐为该作业创建轻量环境
`infercnv_r`，比完整 `epn2_r` 更快：

```bash
source /project/$USER/miniforge3/bin/activate
cd /project/$USER/wslproject/projectmd/hpc
bash 01_prepare_infercnv_env.sh
```

`30_infercnv_gte009.slurm` 会优先使用 `infercnv_r`；如果你已经有完整
`epn2_r`，也会自动 fallback 到 `epn2_r`。

还需要准备 hg38 gene order。建议在登录节点或本地下载 GENCODE GTF 后生成：

```bash
cd /project/$USER/wslproject
mkdir -p output/reference
python projectmd/fig02_infercnv_build_gene_order.py \
    --gtf /project/$USER/reference/gencode.v44.annotation.gtf.gz \
    --out output/reference/gencode_hg38_gene_order.tsv
```

提交 inferCNV：

```bash
cd /project/$USER/wslproject/projectmd/hpc
sbatch 30_infercnv_gte009.slurm
```

监控与拉回结果：

```bash
squeue -u $USER
tail -f /project/$USER/wslproject/logs/infercnv_gte009_<jobid>.out

rsync -avzP \
    zhanglichuan@swarm-xfe.whu.edu.cn:/project/zhanglichuan/wslproject/output/fig2_infercnv_strict/ \
    /home/zlcmc/wslproject/output/fig2_infercnv_strict/
```

如果该作业完成，下一步在本地补 DYNC2H1 gene-level CNV volcano 和 chr11 CNV
heatmap；如果 OOM 或超时，先保留日志，不要直接重跑，优先检查 `MaxRSS` 和
`infercnv_gte009_<jobid>.err`。

## 常见问题

- **作业卡 PD (pending) 很久**: hpxg 队列经常排队, `squeue -u $USER` 看 NODELIST 列, 写 (Resources) 就是在等机器, 写 (Priority) 是在等优先级。
- **OOM (out of memory) killed**: 35k 细胞接近 90 GB 上限。如挂掉, 改用 `10_pyscenic_grn_subsample.slurm` (降采样到 15k) 重跑。
- **超 24 h**: `--num_workers` 拉到 16 一般够。如还不够, 必须降采样或申请 hpib/fat 分区权限。
- **想看运行中的输出**: `tail -f logs/scenic_grn_<jobid>.out`
- **想取消作业**: `scancel <jobid>`
