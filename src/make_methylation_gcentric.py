"""CLI for generating gene-centric methylation matrices using MONTI utilities."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from monti import make_methylation_gcentric


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert probe-level methylation data to a gene-centric matrix."
    )
    parser.add_argument(
        "--methylation",
        required=True,
        type=Path,
        help="Path to the probe-level methylation CSV file produced by the TCGA pipeline.",
    )
    parser.add_argument(
        "--expression",
        required=True,
        type=Path,
        help="Path to the gene expression CSV file aligned to the methylation samples.",
    )
    parser.add_argument(
        "--gene-transcript-map",
        required=True,
        type=Path,
        help="Path to the Ensembl gene-to-transcript mapping table (tab-delimited).",
    )
    parser.add_argument(
        "--promoter-probes",
        required=True,
        type=Path,
        help="Path to the promoter probe annotation table (tab-delimited).",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Destination path for the gene-centric methylation matrix (TSV).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    for path in (args.methylation, args.expression, args.gene_transcript_map, args.promoter_probes):
        if not path.exists():
            raise FileNotFoundError(path)

    df = make_methylation_gcentric(
        str(args.methylation),
        str(args.expression),
        str(args.gene_transcript_map),
        str(args.promoter_probes),
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
