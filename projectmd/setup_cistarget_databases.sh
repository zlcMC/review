#!/usr/bin/env bash
# 下载 pySCENIC 所需的 cisTarget 数据库 (hg19, v9 — 与论文 CellRanger hg19-3.0.0 一致)
# 用法：bash projectmd/setup_cistarget_databases.sh
# 总计 ~1 GB；下载完后跑 conda run -n epn2 python projectmd/fig06_pyscenic_run.py

set -euo pipefail
cd "$(dirname "$0")/.."
DEST="output/fig6/cistarget"
mkdir -p "$DEST"
cd "$DEST"

download() {
    local url="$1" name="$2"
    if [[ -s "$name" ]]; then
        echo "[skip] $name 已存在 ($(du -h "$name" | cut -f1))"
        return
    fi
    echo "[get ] $url"
    wget -c --tries=10 --timeout=60 --read-timeout=60 \
         "$url" -O "$name.part" && mv "$name.part" "$name"
    echo "[ok  ] $name ($(du -h "$name" | cut -f1))"
}

# 1) ranking DB (hg19, 500bp upstream, 7 species, mc9nr) —— ~1 GB
download "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather" \
         "hg19-500bp-upstream-7species.mc9nr.feather"

# 2) motif → TF 注释表
download "https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
         "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

# 3) TF list
download "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt" \
         "allTFs_hg38.txt"

echo
echo '====== 文件清单 ======'
ls -lh

cat <<'EOF'

下载完成。下一步：
  conda run -n epn2 python projectmd/fig06_pyscenic_run.py

断点续传：直接重跑本脚本 (-c)。
代理：HTTPS_PROXY=http://127.0.0.1:7890 bash projectmd/setup_cistarget_databases.sh
EOF
