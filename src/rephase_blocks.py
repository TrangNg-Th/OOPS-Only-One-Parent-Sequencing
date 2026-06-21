#!/usr/bin/env python3

import os
import re
import sys
import argparse
import numpy as np
import pandas as pd


# =============================================================================
# Script to rephase haplotype blocks based on mismatch patterns in parent-child trio data.
# =============================================================================


# =============================================================================
# Argument Parsing
# =============================================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description="Rephase haplotype blocks based on mismatch patterns in parent-child data"
    )

    # Required arguments
    parser.add_argument("--in_block_summary_file", required=True,
                        help="TSV file with summarized mismatches per haplotype block")

    parser.add_argument("--in_ht_file", required=True,
                        help="TSV file with phased SNPs and genotypes per position")

    parser.add_argument("--out_dir", required=True,
                        help="Output directory")


    parser.add_argument("--parent_id", required=True,
                        help="Parent sample ID")

    parser.add_argument("--child_id", required=True,
                        help="Child sample ID")

    parser.add_argument("--window_size", type=int, required=True,
                        help="Window size for analysis")

    # Optional arguments
    parser.add_argument("--threshold", type=int, default=100,
                        help="Threshold value (default: 100)")


    return parser.parse_args()


# =============================================================================
# Main
# =============================================================================
args = parse_args()

in_block_summary_file = args.in_block_summary_file
in_ht_file            = args.in_ht_file
out_dir               = args.out_dir
parent_id             = args.parent_id
child_id              = args.child_id
WINDOW                = args.window_size
threshold             = args.threshold



# =============================================================================
# Setup Output
# =============================================================================
os.makedirs(out_dir, exist_ok=True)

print("=" * 60)
print("Running script to rephase blocks that likely have phase switch errors")
print("=" * 60)

print("Parameters:")
for k, v in vars(args).items():
    print(f"{k:30s}: {v}")

# =============================================================================
# Load block summary data to extract the candidate blocks
# =============================================================================
block_summary_df = pd.read_csv(in_block_summary_file, sep="\t", dtype=str)
ht_df = pd.read_csv(in_ht_file, sep="\t", dtype=str)
ht_df["pos"] = ht_df["pos"].astype(int)

# Results list to store the regions that need to be rephased
regions_to_rephase = []


# Check for blocks with high mismatch counts and rephase them
for block in block_summary_df["PS_child"].unique():
    
    tmp               = block_summary_df[block_summary_df["PS_child"] == block]
    min_n_mismmatches = tmp[["n_mismatches_h0", "n_mismatches_h1"]].astype(int).min(axis=1).iloc[0]
    max_n_mismatches  = tmp[["n_mismatches_h0", "n_mismatches_h1"]].astype(int).max(axis=1).iloc[0]
    block_length      = tmp["block_length"].astype(int).iloc[0]
    
    if int(min_n_mismmatches) >= threshold:
        print(f"Block [{block}] (~{block_length//1000}kbp) - min mismatches : {min_n_mismmatches} >= {threshold} - max mismatches: {max_n_mismatches}")
        print("Proceed to rephase the block.\n")
        
        block_ht        = ht_df[ht_df["PS_child"] == block]
        
        # Find the starting position
        start_of_block  = block_ht["pos"].min()
        end_of_block    = block_ht["pos"].max()
        chrom           = block_ht["chrom"].iloc[0]
        
        # Divide the block into windows, the write out the regions to rephase
        for start in range(start_of_block, end_of_block + 1, WINDOW):
            end = min(start + WINDOW - 1, end_of_block)
            regions_to_rephase.append((chrom, start, end))
            
        
        print(f"Total number of regions to rephase: {len(regions_to_rephase)}")
        
        
        # Write out the regions to rephase to a BED file
        rephase_bed_file = os.path.join(out_dir, f"{parent_id}_{child_id}_rephase_regions.bed")
        
if regions_to_rephase: 
    with open(rephase_bed_file, "w") as f:
        for chrom, start, end in regions_to_rephase:
            f.write(f"{chrom}\t{start}\t{end}\n")   
    print(f"Regions to rephase written to: {rephase_bed_file}")
else:
    print("No regions found that exceed the mismatch threshold. No rephasing needed.")
    

