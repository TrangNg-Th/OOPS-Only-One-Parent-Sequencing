#!/usr/bin/env python3
import pandas as pd
import numpy as np
from collections import defaultdict
import re
import sys
import os

# ------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------
if len(sys.argv) not in (12,13): # one extra arg for passing the python script
    print("Usage: python count_mismatches.py "
          "<input_tsv> <out_dir> <min_dp> <max_dp> "
          "<mom_id> <child_id> <GT_qual> <NV_quantile> "
          "<MM_diff_min> <WINDOW> <RECOUNT> [DNMC_FILE]")

    print()
    print("Arguments:")
    print("  input_tsv     Input TSV file with genotype and phasing information")
    print("  out_dir       Output directory")
    print("  min_dp        Minimum read depth per site")
    print("  max_dp        Maximum read depth per site")
    print("  mom_id        Sample ID of the mother")
    print("  child_id      Sample ID of the child")
    print("  GT_qual       Minimum genotype quality (GQ) per site")
    print(
        "  NV_quantile   Quantile threshold on the number of variants per block; "
        "only blocks with a variant count >= this quantile are considered"
    )
    print(
        "  MM_diff_min   Minimum mismatch fraction between two blocks, defined as "
        "(number of mismatches / number of variants in the block)"
    )
    print(
        "  final_filter (T/F)  Apply a more lose final filtering layer.  "
        "Extract blocks with one mismatch and > 1 or 0 mismatch in the other block"
    )
    print(
        " WINDOW    the window around which the snp was considered. Usually we take all"
    )
    print(
        "recount (T or F) if T (true), the recounting will allow for two haplotype blocks to not have difference in mismatch ( both can have mismatch = 1)"
    )
    sys.exit(1)


in_file = sys.argv[1]
out_dir = sys.argv[2]
min_dp = int(sys.argv[3])
max_dp = int(sys.argv[4])
mom_id = sys.argv[5]
child_id = sys.argv[6]
GT_qual = int(sys.argv[7]) # 30
NV_QUANTILE = float(sys.argv[8]) 
MM_DIFF_MIN = float(sys.argv[9]) #0.3 # at least 30% of differences
WINDOW=int(sys.argv[10])
RECOUNT=sys.argv[11]
DNMC_FILE=None
if RECOUNT == "T":
    if len(sys.argv) != 13:
        print("Error: DNMC_FILE required when RECOUNT == T")
        sys.exit(1)
    DNMC_FILE = sys.argv[12]



# Make sure the window exists
if WINDOW == 0:
    w = ''
else:
    w = str(WINDOW // 1000) + "kb."


os.makedirs(out_dir, exist_ok=True)

## Print out all the parameters for logging
print("--"*20)
print("Input:", in_file)
print("Output dir:", out_dir)
print("DP range:", min_dp, "-", max_dp)
print("Mother ID:", mom_id)
print("Child ID:", child_id)
print("Genotype Quality threshold:", GT_qual)
print("Number of Variants Quantile threshold:", NV_QUANTILE)
print("Mismatch Difference Minimum:", MM_DIFF_MIN)


print("--"*20)


# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
df = pd.read_csv(in_file, sep="\t", dtype=str)

# Add headers
# df.columns = ["chrom", "pos", "PS_mom", "PS_child", "GT_mom", "GT_child", "DP_mom", "DP_child", "GQ_mom", "GQ_child", "None"]


## --------------------------------------------------------------------
# DP filtering
## --------------------------------------------------------------------
df = df[
    (df["DP_mom"].astype(int).between(min_dp, max_dp)) &
    (df["DP_child"].astype(int).between(min_dp, max_dp))
]

# Inlude only alleles where both mom and child have genotype quality >= 30
df = df[(df['GQ_mom'].astype(int) >= GT_qual) & (df['GQ_child'].astype(int) >= GT_qual)]


# ------------------------------------------------------------------
# Genotype sanity filtering
# ------------------------------------------------------------------
biallelic01 = re.compile(r"^[01][\/|][01]$")
df = df[
    df["GT_mom"].str.match(biallelic01) &
    df["GT_child"].str.match(biallelic01)
]

df["pos"] = df["pos"].astype(int)

# ------------------------------------------------------------------
# Per-PS mismatch counting
# ------------------------------------------------------------------
results = []

grouped = df.groupby(["chrom", "PS_child"])
for (chrom, ps), g in grouped:
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
        "block_length": g["pos"].max() - g["pos"].min() + 1,
        "n_variants": g.shape[0],
        "n_mismatches_h0": h0_mm,
        "n_mismatches_h1": h1_mm
    })

summary_df = pd.DataFrame(results)

# Create the ratio between mismatches
min_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].min(axis=1)
max_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].max(axis=1)

summary_df["mismatch_difference"] = np.where(
    max_mm > 0,
    min_mm / max_mm,
    0.0
)



# ------------------------------------------------------------------
# Threshold 1: large PS blocks only
# ------------------------------------------------------------------
n_variants_series = summary_df["n_variants"]
quantiles = n_variants_series.quantile([0, 0.2, 0.25, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99])

print("n_variants quantiles:")
print(quantiles)
print("Writing file *mismatch.tsv")
summary_file = f"{out_dir}/{mom_id}_{child_id}_mismatch.{w}tsv"
summary_df.to_csv(summary_file, sep="\t", index=False)

# Choose cutoff 
nv_cutoff = quantiles.loc[NV_QUANTILE]
print(f"Number of minimum snps considered within a block is {nv_cutoff }")
summary_df = summary_df[summary_df["n_variants"] >= nv_cutoff]


# ------------------------------------------------------------------
# Threshold 2: exactly one mismatch on one haplotype
# ------------------------------------------------------------------
# (this criterium is different when we RECOUNT (using a smaller window))
if RECOUNT == "T":
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
else:
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


# Sanity filter     
summary_df = summary_df[summary_df["block_length"] > 1]

# ------------------------------------------------------------------
# Threshold 3: strong haplotype asymmetry
# ------------------------------------------------------------------
if RECOUNT != "T":
    summary_df = summary_df[summary_df["mismatch_difference"] < MM_DIFF_MIN]


# ------------------------------------------------------------------
# Identify mismatch positions (candidate DNMs)
# ------------------------------------------------------------------
merged = df.merge(summary_df, on=["chrom", "PS_child"], how="inner")

# load the suggested candidates (only use for RECOUNT)
dnmc_set = set()

if RECOUNT == "T":
    dnmc_file = pd.read_csv(DNMC_FILE, sep="\t", dtype=str)

    dnmc_set = set(
        zip(
            dnmc_file["#[1]CHROM"],
            dnmc_file["[2]POS"].astype(int)
        )
    )

# Save results
dnm_records = []

if RECOUNT != "T":
    for (chrom, ps), g in merged.groupby(["chrom", "PS_child"]):
        v = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
        h0, h1 = v[0].values, v[1].values
        mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

        mask_h0 = (mom == h0[:, None]).any(axis=1)
        mask_h1 = (mom == h1[:, None]).any(axis=1)

        if np.sum(~mask_h0) == 1:
            pos = g.loc[~mask_h0, "pos"].iloc[0]
        elif np.sum(~mask_h1) == 1:
            pos = g.loc[~mask_h1, "pos"].iloc[0]
        else:
            continue

        dnm_records.append({
            "chrom": chrom,
            "pos": pos,
            "PS_child": ps
        })
elif RECOUNT == "T":
    
    for (chrom, ps), g in merged.groupby(["chrom", "PS_child"]):
        
        v = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
        h0, h1 = v[0].values, v[1].values
        mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

        mask_h0 = (mom == h0[:, None]).any(axis=1)
        mask_h1 = (mom == h1[:, None]).any(axis=1)

        if (np.sum(~mask_h0) == 1) and (np.sum(~mask_h1) > 1):
            pos = g.loc[~mask_h0, "pos"].iloc[0]

        elif (np.sum(~mask_h1) == 1) and (np.sum(~mask_h0) > 1):
            pos = g.loc[~mask_h1, "pos"].iloc[0]

        elif (np.sum(~mask_h1) == 1) and (np.sum(~mask_h0) == 1):
            pos_h0 = g.loc[~mask_h0, "pos"].iloc[0]
            pos_h1 = g.loc[~mask_h1, "pos"].iloc[0]

            candidates = [
                p for p in (pos_h0, pos_h1)
                if (chrom, p) in dnmc_set
            ]

            if len(candidates) != 1:
                continue

            pos = candidates[0]

        else:
            continue

        if (chrom, pos) not in dnmc_set:
            continue

        dnm_records.append({
            "chrom": chrom,
            "pos": pos,
            "PS_child": ps
        })
    
dnm_df = pd.DataFrame(dnm_records)
dnm_df.to_csv(f"{out_dir}/{mom_id}_{child_id}_dnmc.{w}tsv", sep="\t", index=False)
print("Writing file *dnmc.tsv")
# BED
bed = pd.DataFrame({
    "chrom": dnm_df["chrom"],
    "start": dnm_df["pos"]-1,
    "end": dnm_df["pos"]
}).drop_duplicates()

print("Writing file *dnmc.bed")
bed.to_csv(f"{out_dir}/{mom_id}_{child_id}_dnmc.{w}bed", sep="\t", index=False, header=False)


