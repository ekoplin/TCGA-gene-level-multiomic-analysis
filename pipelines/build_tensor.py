#!/usr/bin/env python3
"""Build a MOPA-compatible tensor from gene-level omics matrices."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import yaml


def read_config(path: Path) -> Dict:
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Matrix not found: {path}")
    df = pd.read_csv(path, sep="\t")
    df = df.set_index(df.columns[0])
    df.index.name = "Gene"
    return df


def harmonise_matrices(
    expr: pd.DataFrame, meth: pd.DataFrame, mir: pd.DataFrame
) -> Tuple[pd.Index, pd.Index, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes = expr.index.intersection(meth.index).intersection(mir.index)
    if genes.empty:
        raise ValueError("No shared genes across matrices. Check preprocessing results.")
    samples = expr.columns.intersection(meth.columns).intersection(mir.columns)
    if samples.empty:
        raise ValueError("No shared samples across matrices. Check preprocessing results.")
    genes = genes.sort_values()
    samples = samples.sort_values()
    expr = expr.loc[genes, samples]
    meth = meth.loc[genes, samples]
    mir = mir.loc[genes, samples]
    return genes, samples, expr, meth, mir


def build_tensor(expr_path: Path, meth_path: Path, mir_path: Path, output_dir: Path) -> None:
    expr = load_matrix(expr_path)
    meth = load_matrix(meth_path)
    mir = load_matrix(mir_path)

    genes, samples, expr, meth, mir = harmonise_matrices(expr, meth, mir)

    tensor = np.stack([expr.values, meth.values, mir.values])
    output_dir.mkdir(parents=True, exist_ok=True)

    tensor_path = output_dir / "tensor.npy"
    np.save(tensor_path, tensor)

    genelist_path = output_dir / "genelist.txt"
    np.savetxt(genelist_path, genes.values, fmt="%s")

    samplelist_path = output_dir / "samplelist.txt"
    np.savetxt(samplelist_path, samples.values, fmt="%s")

    print(f"Saved tensor to {tensor_path}")
    print(f"Saved gene list to {genelist_path}")
    print(f"Saved sample list to {samplelist_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a TCGA multi-omics tensor")
    parser.add_argument(
        "--config",
        default="config/config.yml",
        type=Path,
        help="Path to pipeline configuration file.",
    )
    args = parser.parse_args()

    config = read_config(args.config)
    output_dir = Path(config.get("output_dir", "output"))
    matrices = config.get("matrices", {})

    expr_path = Path(matrices.get("expression", output_dir / "gene_expression_matrix.tsv"))
    meth_path = Path(matrices.get("methylation", output_dir / "gene_methylation_matrix.tsv"))
    mir_path = Path(matrices.get("mirna", output_dir / "gene_miRNA_matrix.tsv"))

    build_tensor(expr_path, meth_path, mir_path, output_dir)


if __name__ == "__main__":
    main()
