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
PS_blocks = pd.read_csv(out_dir + f"/{mom_id}_{child_id}_mismatch.tsv", sep="\t", dtype=str)


## Extract only blocks that passed the filters