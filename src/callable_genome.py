#!/usr/bin/env python3
import pandas as pd
import numpy as np
from collections import defaultdict
import re
import sys
import os


## ================================================================
## Notes about the function
## ================================================================
## Here we extract the PS ids of all the blocks used for calculation of dnm
## We try to estimate the callable genome size



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

## ==================================== Recount ==================================
if RECOUNT == "T":
    if len(sys.argv) != 13:
        print("Error: DNMC_FILE required when RECOUNT == T")
        sys.exit(1)
    DNMC_FILE = sys.argv[12]


## ===================================== Make sure the window exists ====================
if WINDOW == 0:
    w = ''
else:
    w = str(WINDOW // 1000) + "kb."
    


## ===================================== Create directory ================================
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
print("Applying all the filters to all the PS blocks to calculate number of bases")
print("--"*20)


## =================================================================
## Load data
## =================================================================
df = pd.read_csv(in_file, sep="\t", dtype=str)

## Filter by read depth
df = df[
    (df["DP_mom"].astype(int).between(min_dp, max_dp)) &
    (df["DP_child"].astype(int).between(min_dp, max_dp))
]

## Filter by Genotype quality GQ >= GT_qual
df = df[(df['GQ_mom'].astype(int) >= GT_qual) & (df['GQ_child'].astype(int) >= GT_qual)]

## Filter either mono allelic of bi allelic
biallelic01 = re.compile(r"^[01][\/|][01]$")
df = df[df["GT_mom"].str.match(biallelic01) & df["GT_child"].str.match(biallelic01)]

## Set the type to be integer
df["pos"] = df["pos"].astype(int)


## =========================================================================================
## Create a summary dataframe that contains the blocks and their size and the mismatches
## =========================================================================================

def summarize_PS_groups(df, groups):
    # e.g. grouped = df.groupby(["chrom", "PS_child"])
    results = [] ## output list
    grouped = df.groupby(by=groups)
    
    for (chrom, ps), g in grouped:
        v = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
        h0, h1 = v[0].values, v[1].values

        mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

        mask_h0 = (mom == h0[:, None]).any(axis=1)
        mask_h1 = (mom == h1[:, None]).any(axis=1)

        h0_mm = np.sum(~mask_h0)
        h1_mm = np.sum(~mask_h1)
        
        if h0_mm < h1_mm:
            dnmc_hap = "h0"
        else:
            dnmc_hap = "h1"

        results.append({
            "chrom": chrom,
            "PS_child": ps,
            "block_length": g["pos"].max() - g["pos"].min() + 1,
            "n_variants": g.shape[0],
            "n_mismatches_h0": h0_mm,
            "n_mismatches_h1": h1_mm,
            "cand_hapl" : dnmc_hap
        })

    result_df = pd.DataFrame(results)
    return result_df

## ===================== Build PS blocks ==========================
summary_df = summarize_PS_groups(df, ["chrom", "PS_child"])
print(f"BEFORE: Total number of bases {summary_df['block_length'].sum()}")
print(f"BEFORE: Total number of blocks {summary_df['PS_child'].count()}\n")

## Filter by the biggest blocks only
n_variants_series = summary_df["n_variants"]
quantiles = n_variants_series.quantile([0, 0.2, 0.25, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99])
# Choose cutoff 
nv_cutoff = quantiles.loc[NV_QUANTILE]
summary_df = summary_df[summary_df["n_variants"] >= nv_cutoff]
print(f"Only taking blocks that have at least {nv_cutoff} variants")

## Filter the mismatch differences
min_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].min(axis=1)
max_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].max(axis=1)

summary_df["mismatch_difference"] = np.where(max_mm > 0,    min_mm / max_mm, 0.0)
summary_df = summary_df[summary_df["mismatch_difference"] < MM_DIFF_MIN]


## Only take blocks that have at most 1 mismatch
summary_df = summary_df[
        (
            summary_df["n_mismatches_h0"].isin([0, 1]) &
            (summary_df["n_mismatches_h1"] >= 1)
        ) |
        (
            (summary_df["n_mismatches_h0"] >= 1) &
            summary_df["n_mismatches_h1"].isin([0, 1])
        )
]

## Only extract blocks that have at least more than a base
summary_df = summary_df[summary_df["block_length"] > 1]

## Results
print(f"Total number of bases {summary_df['block_length'].sum()}")
print(f"Total number of blocks {summary_df['PS_child'].count()}")


## Write out the PS file
# print("Writing file *mismatch.tsv")
summary_file = f"{out_dir}/{mom_id}_{child_id}_mismatch.{w}tsv"
summary_df.to_csv(summary_file, sep="\t", index=False)

## Write out the original dataframe
ps_file = f"{out_dir}/{mom_id}_{child_id}_ps_filtered.{w}tsv"
df = df[df["PS_child"].isin(summary_df['PS_child'].unique())]


df = pd.merge(left=df, 
              right=summary_df[["chrom", "PS_child", "cand_hapl"]], 
              on=["chrom", "PS_child"], how="left")

## Extract the sites where the mom is 0/1 or 1/0 
## And the child is 0|1
df = df[
    (df["cand_hapl"].isin(["h0"]) & df["GT_child"].isin(["1|0"]) & df["GT_mom"].isin(["1/0", "0/1"]))
    |
     (df["cand_hapl"].isin(["h1"]) & df["GT_child"].isin(["0|1"]) & df["GT_mom"].isin(["1/0", "0/1"]))
]

print(f"Number of sites {df.shape[0]}")


## ====================================================
## Per each PS_child group, sample out one position randomly
## ====================================================
# Set random seed for reproducibility (optional but recommended)
RANDOM_SEED = 42

sampled_df = (
    df.groupby("PS_child", group_keys=False)
      .apply(lambda x: x.sample(1, random_state=RANDOM_SEED) if len(x) > 0 else x)
      .reset_index(drop=True)
)

print(f"Number of sampled blocks {sampled_df.shape[0]}")
bed = pd.DataFrame({
    "chrom": sampled_df["chrom"],
    "start": sampled_df["pos"]-1,
    "end": sampled_df["pos"]
}).drop_duplicates()

print("Writing file *hetc.bed")
bed.to_csv(f"{out_dir}/{mom_id}_{child_id}_hetc.{w}bed", sep="\t", index=False, header=False)
