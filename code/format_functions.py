#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Theo Portlock
Prepare functional microbiome tables for MaAsLin3 analysis.
"""

import pandas as pd
import re
from pathlib import Path


def make_r_safe(names):
    """
    Convert column names into syntactically valid R identifiers.
    Similar to R's make.names().
    """
    clean = []
    seen = {}

    for name in names:
        n = str(name)

        # Replace non-alphanumeric characters with underscore
        n = re.sub(r"[^\w]", "_", n)

        # Prevent names starting with numbers
        if re.match(r"^\d", n):
            n = "X" + n

        # Collapse multiple underscores
        n = re.sub(r"_+", "_", n)

        # Remove leading/trailing underscores
        n = n.strip("_")

        # Ensure uniqueness
        if n in seen:
            seen[n] += 1
            n = f"{n}_{seen[n]}"
        else:
            seen[n] = 0

        clean.append(n)

    return clean


def process_table(infile, outfile, min_prevalence=5):
    """
    Load, clean, and export microbiome functional tables.
    """

    print(f"Processing: {infile}")

    df = pd.read_csv(infile, index_col=0).T

    # Remove zero samples
    df = df.loc[df.sum(axis=1) != 0]

    # Remove zero features
    df = df.loc[:, df.sum(axis=0) != 0]

    # Remove rare features (optional but recommended)
    if min_prevalence:
        df = df.loc[:, (df > 0).sum(axis=0) >= min_prevalence]

    # Make feature names R-safe
    df.columns = make_r_safe(df.columns)

    # Ensure sample ID column name
    df.index.name = "sampleID"

    # Save
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(outfile, sep="\t")

    print(f"Written: {outfile}")
    print(f"Samples: {df.shape[0]}  Features: {df.shape[1]}")
    print()


def main():

    tables = {
        "data/gutcard.csv": "results/gutcard.tsv",
        "data/gutkegg.csv": "results/gutkegg.tsv",
        "data/gutcazy.csv": "results/gutcazy.tsv",
        "data/gutpfam.csv": "results/gutpfam.tsv",
    }

    for infile, outfile in tables.items():
        process_table(infile, outfile)


if __name__ == "__main__":
    main()
