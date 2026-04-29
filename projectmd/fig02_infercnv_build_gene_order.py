#!/usr/bin/env python3
"""Build an inferCNV gene-order table from a GENCODE/Ensembl GTF file."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def parse_attributes(attribute_text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for field in attribute_text.rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        if len(parts) != 2:
            continue
        key, value = parts
        attributes[key] = value.strip().strip('"')
    return attributes


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def build_gene_order(gtf_path: Path, output_path: Path) -> None:
    gene_ranges: dict[str, tuple[str, int, int]] = {}
    allowed_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    allowed_chromosomes |= {str(i) for i in range(1, 23)} | {"X", "Y"}

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chromosome = fields[0]
            if chromosome not in allowed_chromosomes:
                continue
            attributes = parse_attributes(fields[8])
            gene_name = attributes.get("gene_name") or attributes.get("gene_id")
            if not gene_name:
                continue
            start = int(fields[3])
            end = int(fields[4])
            previous = gene_ranges.get(gene_name)
            if previous is None:
                gene_ranges[gene_name] = (chromosome, start, end)
            else:
                prev_chromosome, prev_start, prev_end = previous
                if prev_chromosome == chromosome:
                    gene_ranges[gene_name] = (
                        chromosome,
                        min(prev_start, start),
                        max(prev_end, end),
                    )

    chromosome_rank = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24}
    chromosome_rank |= {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24}
    rows = sorted(
        gene_ranges.items(),
        key=lambda item: (chromosome_rank[item[1][0]], item[1][1], item[1][2], item[0]),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt") as handle:
        for gene_name, (chromosome, start, end) in rows:
            if not chromosome.startswith("chr"):
                chromosome = f"chr{chromosome}"
            handle.write(f"{gene_name}\t{chromosome}\t{start}\t{end}\n")

    print(f"Wrote {len(rows):,} genes to {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, type=Path, help="GENCODE/Ensembl GTF, optionally .gz")
    parser.add_argument("--out", required=True, type=Path, help="inferCNV gene-order TSV")
    args = parser.parse_args()
    build_gene_order(args.gtf, args.out)


if __name__ == "__main__":
    main()