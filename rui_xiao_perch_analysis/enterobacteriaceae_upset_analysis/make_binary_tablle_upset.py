#!/usr/bin/env python3
import pandas as pd

# Load extracted species
df = pd.read_csv("enterobacteriaceae_species.tsv", sep="\t")

# Optional: apply a presence threshold (e.g., >0 reads, or require >0.01%)
# df = df[df["percent"] > 0.01]

df = df[df["reads"] > 0]
df["present"] = 1

# Pivot to binary matrix: rows=species, columns=samples
binary = df.pivot_table(
    index="species",
    columns="sample",
    values="present",
    aggfunc="max",
    fill_value=0,
)

binary.to_csv("enterobacteriaceae_presence_matrix.tsv", sep="\t")
