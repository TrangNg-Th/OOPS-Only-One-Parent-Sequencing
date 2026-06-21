#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import re

# ----------------------------
# Args
# ----------------------------
if len(sys.argv) != 10:
    print("Usage: count_rephase_mismatches.py "
          "<input_tsv> <out_dir> <min_dp> <max_dp> "
          "<parent_id> <child_id> <GT_qual> <nv_quantile> <mm-diff-min>")
    sys.exit(1)

in_file, out_dir, min_dp, max_dp, parent_id, child_id, GT_qual, NV_QUANTILE, MM_DIFF_MIN = sys.argv[1:]
min_dp = int(min_dp)
max_dp = int(max_dp)
GT_qual = int(GT_qual)
NV_QUANTILE = float(NV_QUANTILE)
MM_DIFF_MIN = float(MM_DIFF_MIN)

os.makedirs(out_dir, exist_ok=True)

print("==== Rephase mismatch counting ====")
print("Input:", in_file)

summary_file = f"{out_dir}/{parent_id}_{child_id}_mismatch_rephase.tsv"
tsv_out = f"{out_dir}/{parent_id}_{child_id}_dnmc_rephase.tsv"
bed_out = f"{out_dir}/{parent_id}_{child_id}_dnmc_rephase.bed"

SUMMARY_COLS = ["chrom", "PS_child", "n_variants", "block_length",
                "n_mismatches_h0", "n_mismatches_h1"]


def write_empty_summary_and_exit(msg):
    """Write an empty (header-only) summary file and exit cleanly."""
    print(msg)
    pd.DataFrame(columns=SUMMARY_COLS).to_csv(summary_file, sep="\t", index=False)
    print("Wrote (empty):", summary_file)
    sys.exit(0)


# ----------------------------
# Load
# ----------------------------
df = pd.read_csv(in_file, sep="\t", dtype=str)

if df.empty:
    write_empty_summary_and_exit("Input is empty after load; nothing to count.")

df["pos"] = df["pos"].astype(int)

# ----------------------------
# Basic filtering
# ----------------------------
df = df[
    df["DP_mom"].astype(int).between(min_dp, max_dp) &
    df["DP_child"].astype(int).between(min_dp, max_dp)
]

df = df[
    (df["GQ_mom"].astype(int) >= GT_qual) &
    (df["GQ_child"].astype(int) >= GT_qual)
]

biallelic = re.compile(r"^[01][\/|][01]$")
df = df[
    df["GT_mom"].str.match(biallelic) &
    df["GT_child"].str.match(biallelic)
]

if df.empty:
    write_empty_summary_and_exit("No records remain after depth/quality/biallelic filtering.")

# ----------------------------
# Count mismatches per PS
# ----------------------------
results = []

for (chrom, ps), g in df.groupby(["chrom", "PS_child"]):

    v = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
    h0, h1 = v[0].values, v[1].values

    mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

    mask_h0 = (mom == h0[:, None]).any(axis=1)
    mask_h1 = (mom == h1[:, None]).any(axis=1)

    h0_mm = np.sum(~mask_h0)
    h1_mm = np.sum(~mask_h1)

    results.append({
        "chrom": chrom,
        "PS_child": ps,
        "n_variants": len(g),
        "block_length": g["pos"].max() - g["pos"].min() + 1,
        "n_mismatches_h0": h0_mm,
        "n_mismatches_h1": h1_mm
    })

if len(results) == 0:
    write_empty_summary_and_exit("No PS blocks formed; nothing to count.")

summary_df = pd.DataFrame(results)

# ------------------------------------------------------------------
# Threshold 1: large PS blocks only
# ------------------------------------------------------------------
n_variants_series = summary_df["n_variants"]
quantiles = n_variants_series.quantile([0, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99])

print("n_variants quantiles:")
print(quantiles)
print("Writing file *mismatch.tsv")

# Save mismatch summary (full, before block-level filters)
summary_df.to_csv(summary_file, sep="\t", index=False)

# Choose cutoff
nv_cutoff = quantiles.loc[NV_QUANTILE]
print(f"Number of minimum snps considered within a block is {nv_cutoff }")
summary_df = summary_df[summary_df["n_variants"] >= nv_cutoff]

if summary_df.empty:
    print("No blocks pass the n_variants quantile cutoff. No DNMs found.")
    sys.exit(0)


# ---------------------------
# Create the ratio between mismatches
# ---------------------------
min_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].min(axis=1)
max_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].max(axis=1)

summary_df["mismatch_difference"] = np.where(
    max_mm > 0,
    min_mm / max_mm,
    0.0
)

# ---------------------------------------------------------
# Threshold 2: only one mismatch on one haplotype
# ---------------------------------------------------------
summary_df = summary_df[
    (
        (summary_df["n_mismatches_h0"] == 1) &
        (summary_df["n_mismatches_h1"] >= 1)
    ) |
    (
        (summary_df["n_mismatches_h0"] >= 1) &
        (summary_df["n_mismatches_h1"] == 1)
    )
]


# ----------------------------------------------------
# Threshold 3: haplotype symmetry
# -----------------------------------------------------
summary_df = summary_df[summary_df["mismatch_difference"] < MM_DIFF_MIN]
print(f"Number of PS blocks after filtering by minimum of mismatches difference : {len(summary_df)}")

# Sanity filter
summary_df = summary_df[summary_df["block_length"] > 1]

if summary_df.empty:
    print("No blocks pass the mismatch-structure filters. No DNMs found.")
    sys.exit(0)

# -----------------------------------------------------
# Identify mismatch positions (find candidate DNMs)
# -----------------------------------------------------
merged = df.merge(summary_df, on=["chrom", "PS_child"], how="inner")

# load the suggested candidates (only use for RECOUNT)
dnmc_set = set()

# Save results
dnm_records = []


# ----------------------------
# Identify candidate DNMs
# (simpler rule for rephase)
# ----------------------------
for (chrom, ps), g in merged.groupby(["chrom", "PS_child"]):
    v = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
    h0, h1 = v[0].values, v[1].values
    mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

    mask_h0 = (mom == h0[:, None]).any(axis=1)
    mask_h1 = (mom == h1[:, None]).any(axis=1)

    n0 = np.sum(~mask_h0)
    n1 = np.sum(~mask_h1)

    # Strict asymmetry, matching count_mismatches.py RECOUNT:
    # the candidate (DNM-bearing) haplotype has exactly ONE mismatch, and the
    # OTHER (inherited) haplotype has strictly MORE than one. This excludes the
    # (1, 0) clean-het case, where one haplotype has a single stray mismatch and
    # the other is perfectly parent-consistent -- that is not a DNM signal.
    if (n0 == 1) and (n1 > 1):
        pos = g.loc[~mask_h0, "pos"].iloc[0]
    elif (n1 == 1) and (n0 > 1):
        pos = g.loc[~mask_h1, "pos"].iloc[0]
    else:
        continue

    dnm_records.append({
        "chrom": chrom,
        "pos": pos,
        "PS_child": ps
    })

if len(dnm_records) == 0:
    print("No DNMs found")
    sys.exit(0)

dnm_df = pd.DataFrame(dnm_records)
print(f"Number of candidate DNMs identified (before LR validation): {len(dnm_df)}")

dnm_df.to_csv(tsv_out, sep="\t", index=False)

bed = pd.DataFrame({
    "chrom": dnm_df["chrom"],
    "start": dnm_df["pos"] - 1,
    "end": dnm_df["pos"]
}).drop_duplicates()

bed.to_csv(bed_out, sep="\t", index=False, header=False)

print("Wrote:", tsv_out)
print("Wrote:", bed_out)