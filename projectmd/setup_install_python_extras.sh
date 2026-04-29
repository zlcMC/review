#!/usr/bin/env bash
# 复现 fimmu-13-903246 所需的 Python 工具（在 epn2 环境内安装）
# 用法: bash projectmd/setup_install_python_extras.sh
set -euo pipefail

PY=/home/zlcmc/miniconda3/envs/epn2/bin/python
PIP="$PY -m pip"

# pip 26.0.1 有 from_json bug，固定到 24.x
$PIP install --quiet "pip<25"

# velocyto 编译时需要 numpy/cython 已存在
$PIP install --quiet "numpy<2" cython

# velocyto: 必须 --no-build-isolation 才能找到上面装的 numpy
$PIP install --quiet --no-build-isolation velocyto

# scVelo + pySCENIC
$PIP install --quiet scvelo pyscenic

# 可选: cellrank (scVelo 下游)
$PIP install --quiet cellrank || true

echo "Done. Verify:"
$PY - <<'PYEOF'
import importlib, sys
for pkg in ["scvelo", "velocyto", "pyscenic"]:
    try:
        m = importlib.import_module(pkg)
        v = getattr(m, "__version__", "?")
        print(f"  {pkg:10s} {v}")
    except Exception as e:
        print(f"  {pkg:10s} FAIL: {e}", file=sys.stderr)
PYEOF
