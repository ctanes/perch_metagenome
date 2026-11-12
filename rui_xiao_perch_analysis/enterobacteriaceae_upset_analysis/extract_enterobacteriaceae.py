#!/usr/bin/env python3
"""
extract_enterobacteriaceae.py
Extract species (rank S) under Enterobacteriaceae family from Kraken2 reports.
Output: sample  species  taxid  reads  percent
"""

import glob
import sys
from pathlib import Path


def extract_species(report_path):
    sample = Path(report_path).stem.replace("_k2", "").replace("_k2.txt", "")
    entries = []
    in_family = False
    family_depth = None

    with open(report_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 6:
                continue

            percent, num_clade, num_direct, rank, taxid = parts[:5]
            name_field = parts[5]

            # Depth from leading spaces (Kraken2 uses 2 spaces per level, but we just count raw spaces)
            stripped = name_field.lstrip(" ")
            depth = len(name_field) - len(stripped)
            name = stripped

            # Enterobacteriaceae family line
            if rank == "F" and name == "Enterobacteriaceae":
                in_family = True
                family_depth = depth
                continue

            # Exited the family subtree (same or shallower depth in DFS)
            if in_family and depth <= family_depth:
                in_family = False
                family_depth = None
                continue

            # Species under Enterobacteriaceae
            if in_family and rank == "S":
                entries.append({
                    "sample": sample,
                    "species": name,
                    "taxid": taxid,
                    "reads": int(num_clade),
                    "percent": float(percent),
                })

    return entries


if __name__ == "__main__":
    reports = glob.glob("*_k2.txt")
    if not reports:
        print("No reports found matching *_k2.txt", file=sys.stderr)
        sys.exit(1)

    print("sample\tspecies\ttaxid\treads\tpercent")
    for report in reports:
        for e in extract_species(report):
            print(f"{e['sample']}\t{e['species']}\t{e['taxid']}\t{e['reads']}\t{e['percent']:.4f}")