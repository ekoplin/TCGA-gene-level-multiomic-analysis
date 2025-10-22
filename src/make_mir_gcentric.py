"""CLI for generating gene-centric miRNA matrices using MONTI utilities."""

from __future__ import annotations

import argparse
from pathlib import Path

from monti import make_mir_gcentric


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert mature miRNA expression to a gene-centric matrix using MONTI mappings."
    )
    parser.add_argument(
        "--mirna",
        required=True,
        type=Path,
        help="Path to the miRNA expression CSV file produced by the TCGA pipeline.",
    )
    parser.add_argument(
        "--expression",
        required=True,
        type=Path,
        help="Path to the gene expression CSV file aligned to the miRNA samples.",
    )
    parser.add_argument(
        "--targets",
        required=True,
        type=Path,
        help="Path to the miRNA target mapping table (tab-delimited).",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Destination path for the gene-centric miRNA matrix (TSV).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    for path in (args.mirna, args.expression, args.targets):
        if not path.exists():
            raise FileNotFoundError(path)

    df = make_mir_gcentric(
        str(args.mirna),
        str(args.expression),
        str(args.targets),
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
