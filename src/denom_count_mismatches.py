#!/usr/bin/env python3
import pandas as pd
import numpy as np
from collections import defaultdict
import re
import sys
import os


## In this function, we will extract the PS blocks IDs that have 0 or 1 mismatches.
## Initially, I'll apply the same filters to extract the blocks,
## Then, I'll extract the block directly from HTblocks/*ps.tsv
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


print("--"*20)


# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
df = pd.read_csv(in_file, sep="\t", dtype=str)

## --------------------------------------------------------------------
## STEP 0 : FILTER TO CALCULATE MISMATCH DIFFERENCES TO OBTAIN BLOCKS THAT HAVE
##                                  0 OR 1 MISMATCHES
## --------------------------------------------------------------------
## DP filtering 
df = df[
    (df["DP_mom"].astype(int).between(min_dp, max_dp)) &
    (df["DP_child"].astype(int).between(min_dp, max_dp))
]

## GQ >= GT_qual
df = df[(df['GQ_mom'].astype(int) >= GT_qual) & (df['GQ_child'].astype(int) >= GT_qual)]


# ------------------------------------------------------------------
## Biallelic sites filter
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


def summarize_PS_groups(df, groups):
    # grouped = df.groupby(["chrom", "PS_child"])
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
print(f"Total number of bases {summary_df['block_length'].sum()}")
print(f"Total number of blocks {summary_df['PS_child'].count()}")

## ========================= Create the ratio between mismatches
min_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].min(axis=1)
max_mm = summary_df[["n_mismatches_h0", "n_mismatches_h1"]].max(axis=1)

summary_df["mismatch_difference"] = np.where(
    max_mm > 0,
    min_mm / max_mm,
    0.0
)
n_variants_series = summary_df["n_variants"]
quantiles = n_variants_series.quantile([0, 0.2, 0.25, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99])
# Choose cutoff 
nv_cutoff = quantiles.loc[NV_QUANTILE]
summary_df = summary_df[summary_df["n_variants"] >= nv_cutoff]
summary_df = summary_df[summary_df["mismatch_difference"] < MM_DIFF_MIN]
summary_df = summary_df[summary_df["block_length"] > 1]
# ------------------------------------------------------------------
# Threshold 2: exactly one mismatch on one haplotype ==> APPLY
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
            summary_df["n_mismatches_h0"].isin([0, 1]) &
            (summary_df["n_mismatches_h1"] >= 1)
        ) |
        (
            (summary_df["n_mismatches_h0"] >= 1) &
            summary_df["n_mismatches_h1"].isin([0, 1])
        )
]
    
print(f"The genome size of considered {summary_df['block_length'].sum()}")
print(f"The number of blocks considered {summary_df['PS_child'].count()}")
print()
## =================================================================================
## Part 2: Extract the PS blocks from the child
##         and save it
## ==================================================================================
PS_blocks = summary_df['PS_child'].unique()

## Load again the original file
df = pd.read_csv(in_file, sep="\t", dtype=str)

## Only retain biallelic sites
biallelic01 = re.compile(r"^[01][\/|][01]$")
df = df[df["GT_mom"].str.match(biallelic01) & df["GT_child"].str.match(biallelic01)]
df["pos"] = df["pos"].astype(int)  ## set POS to be integer

## Extract blocks that are actually from the candidate block
df = df[df['PS_child'].isin(PS_blocks)]

## Also add the columns whihc haplotype to consider
df = pd.merge(left=df, right=summary_df[["chrom", "PS_child", "cand_hapl"]], 
              how="left", on=["chrom", "PS_child"]
              )

# ## ==============  Now only extract the number of hetero in child, homo in mom
# ## ============== Since we only considered the mismatches as the hetero variant, we must extract only
# ## ============== positions where the snps are on the min haplotype
## =================================================================
results = []
grouped = df.groupby(["chrom", "PS_child"])
not_filtered_candidates = 0

for (chrom, ps), g in grouped:
    ## Filters mom homo, child hetero
    g = g[g["GT_mom"].isin(["0/0"]) & g["GT_child"].isin(["0|1", "1|0"])]
    
    if len(g) > 0 :
 
        v = g["GT_child"].str.split(r"|", expand=True).astype(int)
        h0, h1 = v[0].values, v[1].values
        
        
        # # count the number of hetero that passes the filters
        if (g["cand_hapl"] == "h0").all():
            not_filtered_candidates += np.count_nonzero(h0 == 1)
            # print(f"Number of candidates {np.count_nonzero(h0 == 1)}")
        elif (g["cand_hapl"] == "h1").all():
            not_filtered_candidates += np.count_nonzero(h0 == 1)
            # print(f"Number of candidates {np.count_nonzero(h1 == 1)}")
            
print(f"Number of HET bases  for the child, HOMO REF for mum is: {not_filtered_candidates}")
## ============================================================================
# ## Add filterings and see how many have left
## ============================================================================
df = df[
    (df["DP_mom"].astype(int).between(min_dp, max_dp)) &
    (df["DP_child"].astype(int).between(min_dp, max_dp))
]

# # Inlude only alleles where both mom and child have genotype quality >= 30
df = df[(df['GQ_mom'].astype(int) >= GT_qual) & (df['GQ_child'].astype(int) >= GT_qual)]

# # # ------------------------------------------------------------------
# # # Per-PS mismatch counting
# # # ------------------------------------------------------------------
results2 = []

grouped = df.groupby(["chrom", "PS_child"])
for (chrom, ps), g in grouped:
    v = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
    h0, h1 = v[0].values, v[1].values

    mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

    mask_h0 = (mom == h0[:, None]).any(axis=1)
    mask_h1 = (mom == h1[:, None]).any(axis=1)

    h0_mm = np.sum(~mask_h0)
    h1_mm = np.sum(~mask_h1)
    
    if h0_mm < h1_mm : 
        smaller_block = "h0"
    else:
        smaller_block = "h1"

    results2.append({
        "chrom": chrom,
        "PS_child": ps,
        "block_length": g["pos"].max() - g["pos"].min() + 1,
        "n_variants": g.shape[0],
        "n_mismatches_h0": h0_mm,
        "n_mismatches_h1": h1_mm,
        "smaller_block" : smaller_block 
        
    })

summary_df2 = pd.DataFrame(results2)


# # Create the ratio between mismatches
min_mm = summary_df2[["n_mismatches_h0", "n_mismatches_h1"]].min(axis=1)
max_mm = summary_df2[["n_mismatches_h0", "n_mismatches_h1"]].max(axis=1)

summary_df2["mismatch_difference"] = np.where(
    max_mm > 0,
    min_mm / max_mm,
    0.0
)

summary_df2 = summary_df2[
        (
            summary_df2["n_mismatches_h0"].isin([0, 1]) &
            (summary_df2["n_mismatches_h1"] >= 1)
        ) |
        (
            (summary_df2["n_mismatches_h0"] >= 1) &
            summary_df2["n_mismatches_h1"].isin([0, 1])
        )
]


n_variants_series = summary_df2["n_variants"]
quantiles = n_variants_series.quantile([0, 0.2, 0.25, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99])
# Choose cutoff 
nv_cutoff = quantiles.loc[NV_QUANTILE]
summary_df2 = summary_df2[summary_df2["n_variants"] >= nv_cutoff]
summary_df2 = summary_df2[summary_df2["mismatch_difference"] < MM_DIFF_MIN]
summary_df2 = summary_df2[summary_df2["block_length"] > 1]

## ======================== AFTER ALL FILTERS, HOW MANY HET SITES left ? ====
PS_Groups = summary_df2["PS_child"].unique()
results3 = []
grouped3 = df.groupby(["chrom", "PS_child"])
not_filtered_candidates = 0
for (chrom, ps), g in grouped:
    ## Filters mom homo, child hetero
    g = g[g["GT_mom"].isin(["0/0"]) & g["GT_child"].isin(["0|1", "1|0"])]
    
    if len(g) > 0 :
 
        v = g["GT_child"].str.split(r"|", expand=True).astype(int)
        h0, h1 = v[0].values, v[1].values
        
        # # count the number of hetero that passes the filters
        if (g["cand_hapl"] == "h0").all():
            not_filtered_candidates += np.count_nonzero(h0 == 1)
            # print(f"Number of candidates {np.count_nonzero(h0 == 1)}")
        elif (g["cand_hapl"] == "h1").all():
            not_filtered_candidates += np.count_nonzero(h0 == 1)
            # print(f"Number of candidates {np.count_nonzero(h1 == 1)}")
        
print(f"Number of HET bases  for the child, HOMO REF for mum AFTER FILT: {not_filtered_candidates}")


## ============================= Write tsv and bed file to continue the filters ==========

# Save results
df = df[df["PS_child"].isin(PS_Groups)]
df = df[df["GT_mom"].isin(['0/0']) & df["GT_child"].isin(['0|1', '1|0'])]
## Filter for specific rows as well
filtered_df = df[((df["cand_hapl"] == "h0") & (df["GT_child"] == "1|0")) | 
                 ((df["cand_hapl"] == "h1") & (df["GT_child"] == "0|1"))]

## Filter out the columns
filtered_df = filtered_df[["chrom", "pos"]]
filtered_df["start"] = filtered_df["pos"] - 1
filtered_df["end"] = filtered_df["pos"]
filtered_df.drop(["pos"], inplace=True, axis=1)
filtered_df = filtered_df.drop_duplicates()
print("Writing file *hetc.bed")
print(f"Number of HET site passed {filtered_df.shape[0]}")
filtered_df.to_csv(f"{out_dir}/{mom_id}_{child_id}_hetc.{w}bed", sep="\t", index=False, header=False)

