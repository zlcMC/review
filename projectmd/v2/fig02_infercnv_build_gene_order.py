#!/usr/bin/env python3
"""Build an inferCNV gene-order table from a GENCODE/Ensembl GTF file (v2).

v2 改动：仅薄包装。
- 接受可选 --out；缺省写入 output/v2/reference/gene_order.txt
- --gtf 必填（数据通常 >100 MB，不在仓库内）
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from helpers.helpers_common import v2_file


ALLOWED = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
ALLOWED |= {str(i) for i in range(1, 23)} | {"X", "Y"}
RANK = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24}
RANK |= {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24}


def parse_attributes(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in text.rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        if len(parts) == 2:
            out[parts[0]] = parts[1].strip().strip('"')
    return out


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def build_gene_order(gtf: Path, out: Path) -> None:
    ranges: dict[str, tuple[str, int, int]] = {}
    with open_text(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene" or f[0] not in ALLOWED:
                continue
            attrs = parse_attributes(f[8])
            name = attrs.get("gene_name") or attrs.get("gene_id")
            if not name:
                continue
            chrom, start, end = f[0], int(f[3]), int(f[4])
            prev = ranges.get(name)
            if prev is None or prev[0] != chrom:
                ranges[name] = (chrom, start, end)
            else:
                ranges[name] = (chrom, min(prev[1], start), max(prev[2], end))

    rows = sorted(ranges.items(), key=lambda kv: (RANK[kv[1][0]], kv[1][1], kv[1][2], kv[0]))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wt") as fh:
        for name, (chrom, start, end) in rows:
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            fh.write(f"{name}\t{chrom}\t{start}\t{end}\n")
    print(f"Wrote {len(rows):,} genes to {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True, type=Path, help="GENCODE/Ensembl GTF, optionally .gz")
    ap.add_argument("--out", type=Path, default=Path(v2_file("reference", "gene_order.txt")),
                    help="default: output/v2/reference/gene_order.txt")
    a = ap.parse_args()
    build_gene_order(a.gtf, a.out)


if __name__ == "__main__":
    main()
