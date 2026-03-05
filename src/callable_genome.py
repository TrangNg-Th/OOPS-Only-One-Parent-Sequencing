#!/usr/bin/env python3

import os
import re
import sys
import numpy as np
import pandas as pd


# =============================================================================
# Callable Genome Estimation Script
# =============================================================================
# This script:
#   1. Filters phased SNPs by depth / quality
#   2. Builds haplotype blocks (PS groups)
#   3. Filters blocks by size + mismatch structure
#   4. Samples one SNP per qualified block
#   5. Outputs callable genome statistics + BED file
# =============================================================================


# =============================================================================
# Argument Parsing
# =============================================================================
if len(sys.argv) not in (12, 13):
    print("Usage: python count_mismatches.py "
          "<input_tsv> <out_dir> <min_dp> <max_dp> "
          "<mom_id> <child_id> <GT_qual> <NV_quantile> "
          "<MM_diff_min> <WINDOW> <RECOUNT> [DNMC_FILE]")
    sys.exit(1)

in_file      = sys.argv[1]
out_dir      = sys.argv[2]
min_dp       = int(sys.argv[3])
max_dp       = int(sys.argv[4])
mom_id       = sys.argv[5]
child_id     = sys.argv[6]
GT_qual      = int(sys.argv[7])
NV_QUANTILE  = float(sys.argv[8])
MM_DIFF_MIN  = float(sys.argv[9])
WINDOW       = int(sys.argv[10])
RECOUNT      = sys.argv[11]
DNMC_FILE    = None

if RECOUNT == "T":
    if len(sys.argv) != 13:
        print("Error: DNMC_FILE required when RECOUNT == T")
        sys.exit(1)
    DNMC_FILE = sys.argv[12]


# =============================================================================
# Window Label
# =============================================================================
if WINDOW == 0:
    w = ""
else:
    w = f"{WINDOW // 1000}kb."


# =============================================================================
# Setup Output
# =============================================================================
os.makedirs(out_dir, exist_ok=True)

print("=" * 60)
print("Calculating callable genome size")
print("=" * 60)
print("Input file:", in_file)
print("Output dir:", out_dir)
print("Depth range:", min_dp, "-", max_dp)
print("Mother ID:", mom_id)
print("Child ID:", child_id)
print("Genotype Quality:", GT_qual)
print("Variant Quantile:", NV_QUANTILE)
print("Mismatch Difference:", MM_DIFF_MIN)
print("=" * 60)


# =============================================================================
# Load and Initial Filtering
# =============================================================================
df = pd.read_csv(in_file, sep="\t", dtype=str)

# Depth filter
df = df[
    (df["DP_mom"].astype(int).between(min_dp, max_dp)) &
    (df["DP_child"].astype(int).between(min_dp, max_dp))
]

# Genotype quality filter
df = df[
    (df["GQ_mom"].astype(int) >= GT_qual) &
    (df["GQ_child"].astype(int) >= GT_qual)
]

# Keep only biallelic 0/1 genotypes
biallelic01 = re.compile(r"^[01][\/|][01]$")
df = df[
    df["GT_mom"].str.match(biallelic01) &
    df["GT_child"].str.match(biallelic01)
]

df["pos"] = df["pos"].astype(int)


# =============================================================================
# Build PS Block Summary
# =============================================================================
def summarize_PS_groups(df, groups):

    results = []
    grouped = df.groupby(groups)

    for (chrom, ps), g in grouped:

        child_split = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
        h0 = child_split[0].values
        h1 = child_split[1].values

        mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

        mask_h0 = (mom == h0[:, None]).any(axis=1)
        mask_h1 = (mom == h1[:, None]).any(axis=1)

        h0_mm = np.sum(~mask_h0)
        h1_mm = np.sum(~mask_h1)

        dnmc_hap = "h0" if h0_mm < h1_mm else "h1"

        results.append({
            "chrom": chrom,
            "PS_child": ps,
            "block_length": g["pos"].max() - g["pos"].min() + 1,
            "n_variants": g.shape[0],
            "n_mismatches_h0": h0_mm,
            "n_mismatches_h1": h1_mm,
            "cand_hapl": dnmc_hap
        })

    return pd.DataFrame(results)


summary_df = summarize_PS_groups(df, ["chrom", "PS_child"])


# =============================================================================
# Block-Level Filtering
# =============================================================================

# Quantile filter on number of variants
quantiles = summary_df["n_variants"].quantile(
    [0, 0.2, 0.25, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]
)
nv_cutoff = quantiles.loc[NV_QUANTILE]
summary_df = summary_df[summary_df["n_variants"] >= nv_cutoff]

# Mismatch fraction filter
min_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].min(axis=1)
max_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].max(axis=1)

summary_df["mismatch_difference"] = np.where(
    max_mm > 0, min_mm / max_mm, 0.0
)

summary_df = summary_df[
    summary_df["mismatch_difference"] < MM_DIFF_MIN
]

# Keep blocks with <=1 mismatch in one haplotype
summary_df = summary_df[
    (
        summary_df["n_mismatches_h0"].isin([0, 1]) &
        (summary_df["n_mismatches_h1"] >= 1)
    )
    |
    (
        (summary_df["n_mismatches_h0"] >= 1) &
        summary_df["n_mismatches_h1"].isin([0, 1])
    )
]

# Remove trivial blocks
summary_df = summary_df[summary_df["block_length"] > 1]

# Sanity check
# print("Total number of bases considered:",
#       summary_df["block_length"].sum())
# print("Total number of blocks:",
#       summary_df["PS_child"].count())


# ============================================================
# Write callable genome summary file
# ============================================================
callable_file = os.path.join(out_dir, "callable_genome.txt")
total_bases = summary_df["block_length"].sum()

with open(callable_file, "w") as f:
    f.write(f"total_callable_bases\t{total_bases}\n")


# =============================================================================
# Write Block Summary
# =============================================================================
summary_file = f"{out_dir}/{mom_id}_{child_id}_mismatch.{w}tsv"
summary_df.to_csv(summary_file, sep="\t", index=False)


# =============================================================================
# Filter Original DF by Qualified Blocks
# =============================================================================
df = df[df["PS_child"].isin(summary_df["PS_child"].unique())]

df = pd.merge(
    df,
    summary_df[["chrom", "PS_child", "cand_hapl"]],
    on=["chrom", "PS_child"],
    how="left"
)

# Extract candidate heterozygous sites
df = df[
    (
        (df["cand_hapl"] == "h0") &
        (df["GT_child"] == "1|0") &
        (df["GT_mom"].isin(["1/0", "0/1"]))
    )
    |
    (
        (df["cand_hapl"] == "h1") &
        (df["GT_child"] == "0|1") &
        (df["GT_mom"].isin(["1/0", "0/1"]))
    )
]


# =============================================================================
# Sample One SNP Per PS Block
# =============================================================================
RANDOM_SEED = 42

sampled_df = (
    df.groupby("PS_child", group_keys=False)
      .apply(lambda x: x.sample(1, random_state=RANDOM_SEED)
             if len(x) > 0 else x)
      .reset_index(drop=True)
)

print("Number of snps sampled:", sampled_df.shape[0])

total_sampled = sampled_df.shape[0]
with open(callable_file, "a") as f:
    f.write(f"total_sampled_snps\t{total_sampled}\n")



# =============================================================================
# Write BED Output
# =============================================================================
bed = pd.DataFrame({
    "chrom": sampled_df["chrom"],
    "start": sampled_df["pos"] - 1,
    "end": sampled_df["pos"]
}).drop_duplicates()

bed_file = f"{out_dir}/{mom_id}_{child_id}_hetc.{w}bed"
bed.to_csv(bed_file, sep="\t", index=False, header=False)

print("Qualified sampled snps written to bed file:", bed_file)

