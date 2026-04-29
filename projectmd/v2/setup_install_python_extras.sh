#!/usr/bin/env bash
# v2: 复现 fimmu-13-903246 所需的 Python 工具（在 epn2 环境内安装）
# 用法: bash projectmd/v2/setup_install_python_extras.sh
#
# v2 changes:
#   * Idempotency check: 跳过已安装且版本匹配的包，避免重复 pip install
#   * 解除 cellrank 的 `|| true` 静默吞错，改为带提示但不中断
#   * `set -euo pipefail` 已经存在，本脚本沿用

set -euo pipefail

PY=/home/zlcmc/miniconda3/envs/epn2/bin/python
PIP=("$PY" -m pip)

# 工具：检查包是否已安装（任何版本）
have() { "$PY" -c "import importlib, sys; importlib.import_module('$1')" 2>/dev/null; }

# 工具：检查 pip 版本是否符合表达式（这里只判断 < 25）
pip_lt_25() {
    "$PY" - <<'PYEOF' 2>/dev/null
import sys
try:
    from importlib.metadata import version
except Exception:
    from pkg_resources import get_distribution
    v = get_distribution('pip').version
else:
    v = version('pip')
parts = [int(x) for x in v.split('.')[:2] if x.isdigit()]
sys.exit(0 if parts and parts[0] < 25 else 1)
PYEOF
}

# 1) pip 26.0.1 有 from_json bug，固定到 24.x；只在不满足时升降级
if pip_lt_25; then
    echo "[skip] pip already < 25"
else
    "${PIP[@]}" install --quiet "pip<25"
fi

# 2) numpy + cython（velocyto 编译依赖）
if ! have numpy; then "${PIP[@]}" install --quiet "numpy<2"; else echo "[skip] numpy"; fi
if ! have Cython; then "${PIP[@]}" install --quiet cython; else echo "[skip] cython"; fi

# 3) velocyto 必须 --no-build-isolation 才能找到上面装的 numpy
if ! have velocyto; then
    "${PIP[@]}" install --quiet --no-build-isolation velocyto
else
    echo "[skip] velocyto"
fi

# 4) scvelo + pyscenic
for pkg in scvelo pyscenic; do
    if ! have "$pkg"; then "${PIP[@]}" install --quiet "$pkg"; else echo "[skip] $pkg"; fi
done

# 5) 可选 cellrank：失败时给出明确提示，不中断
if ! have cellrank; then
    if ! "${PIP[@]}" install --quiet cellrank; then
        echo "[warn] cellrank install failed (optional, continuing)"
    fi
else
    echo "[skip] cellrank"
fi

echo
echo "Done. Verify:"
"$PY" - <<'PYEOF'
import importlib, sys
for pkg in ["scvelo", "velocyto", "pyscenic", "cellrank"]:
    try:
        m = importlib.import_module(pkg)
        v = getattr(m, "__version__", "?")
        print(f"  {pkg:10s} {v}")
    except Exception as e:
        print(f"  {pkg:10s} MISSING: {e}", file=sys.stderr)
PYEOF
