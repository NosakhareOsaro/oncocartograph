"""GSE96058 gene expression matrix parsing.

The real expression file (confirmed by downloading and inspecting it
before this module was written) is a single ~564MB gzipped CSV: genes as
rows (indexed by **gene symbol**, not Ensembl ID -- unlike this
project's TCGA-derived RNA-seq view), samples as columns (matching the
clinical ``title`` field, e.g. ``"F1"``), values already
log2(FPKM + 0.1)-transformed.

Loading all ~30,865 genes x 3,409 samples into memory is unnecessary when
only a small set of candidate genes is needed for validation -- this
module streams the decompressed file row by row and keeps only rows
whose gene symbol is in the requested set, rather than loading the full
matrix.
"""

from __future__ import annotations

import csv
import gzip
from collections.abc import Iterable
from pathlib import Path

import pandas as pd


def read_selected_gene_expression(path: Path, gene_symbols: Iterable[str]) -> pd.DataFrame:
    """Stream-read only the requested genes from the full GSE96058 expression matrix.

    Args:
        path: Path to the gzipped expression CSV (e.g.
            ``GSE96058_gene_expression_..._transformed.csv.gz``).
        gene_symbols: Gene symbols to extract (rows not in this set are
            skipped without being fully parsed into memory).

    Returns:
        A DataFrame indexed by gene symbol, one column per sample
        (matching the file's header, e.g. ``"F1"``), containing only the
        requested genes that were actually found in the file. Requested
        genes absent from the file are simply not present in the result.
    """
    wanted = set(gene_symbols)
    rows: dict[str, list[float]] = {}
    sample_names: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)
        header = next(reader)
        sample_names = header[1:]
        for row in reader:
            gene = row[0]
            if gene in wanted:
                rows[gene] = [float(x) for x in row[1:]]

    return pd.DataFrame.from_dict(rows, orient="index", columns=sample_names)
