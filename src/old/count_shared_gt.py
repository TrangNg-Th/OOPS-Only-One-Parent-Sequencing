import pandas as pd
import numpy as np
from collections import defaultdict
import re
import os
import sys

# Set working directory to the script's directory
os.chdir(os.path.dirname(os.path.abspath(__file__))) 

# Parse command line arguments
in_file = sys.argv[1]       # input tsv file. Expected : mom_child_PS_fixed.tsv
out_dir = sys.argv[2]       # output directory
min_readdepth = sys.argv[3]
max_readdepth = sys.argv[4]

# Print input parameters
print("Input file:", in_file)
print("Output directory:", out_dir)
print("Min read depth:", min_readdepth)
print("Max read depth:", max_readdepth)

# The paper reported avg coverage for illumina variant call was 30x
min_readdepth = int(min_readdepth)
max_readdepth = int(max_readdepth)

# regex: exactly two alleles (0 or 1) separated by / or |
biallelic01 = re.compile(r'^[01][\/|][01]$')
n_sites = defaultdict(int)
full_data = []

# Load the data
df = pd.read_csv(in_file, sep="\t", dtype=str, header=0)

# Include only alleles where both mom and child have genotype quality 
df = df[((df['DP_mom'].astype(int) >= min_readdepth) & (df['DP_mom'].astype(int) <= max_readdepth)) \
    & ((df['DP_child'].astype(int) >= min_readdepth) & (df['DP_child'].astype(int) <= max_readdepth))]

# Inlude only alleles where both mom and child have genotype quality >= 30
# df = df[(df['GQ_mom'].astype(int) >= 30) & (df['GQ_child'].astype(int) >= 30)]

# Include only homozygous reference or heterozygous genotypes for mom
df = df[df['GT_mom'].isin(['0/0', '1/1'])]

# Write data
df.to_csv(out_dir + "/mom_child_ps_fixed_dp_filtered.tsv", sep="\t", index=False)

# Skip any genotype that is not strictly 0/1 diploid (e.g. 2/2, 0/2, 1/2, 0/0/1, etc.)
df = df[df['GT_mom'].str.match(biallelic01) & df['GT_child'].str.match(biallelic01)]

# Make sure the column "pos" is integer
df['pos'] = df['pos'].astype(int)

# Group by chromosome and phase set
grouped = df.groupby(["chrom", "PS_child"])
for (chrom, ps_child), group in grouped:
    block_length = group["pos"].max() - group["pos"].min() + 1
    n_variants = group.shape[0]
    v = group['GT_child'].str.split(r'[/|]', expand=True).astype(int)
    h0 = v[0].values # haplotype 1
    h1 = v[1].values # haplotype 2
    
    # Turn mom genotypes into a 2D array
    mom_genotypes = group['GT_mom'].str.split(r'[/|]', expand=True).astype(int)
    pairs = mom_genotypes.values  # shape (n, 2)
    
    # Check for matches
    mask_h0 = (pairs == h0[:, None]).any(axis=1)
    mask_h1 = (pairs == h1[:, None]).any(axis=1)
    h0_mismatches = np.sum(1 - mask_h0)
    h1_mismatches = np.sum(1 - mask_h1)
    
    # Save the mismatches counts and total length
    n_sites[(chrom, ps_child)] = (h0_mismatches, h1_mismatches, n_variants, block_length)
    

out_file = out_dir + "/per_chromosome_shared_allele_counts_by_PS.tsv"
# Output results
# print("Writing results to file:", out_file )
# print("Final data shape (to be written):", len(n_sites) )
with open(out_file, "w") as out:
    out.write("chrom\tPS_child\tblock_length\tn_variants\tn_mismatches_h0\tn_mismatches_h1\tinherited_haplotype\tcertainty_\%\n")
    # out.write("chrom\tPS_child\tn_sites\tn_mismatches_h0\tn_mismatches_h1\tpercent_mismatched_h0\tpercent_mismatched_h1\tinherited_haplotype\tcertainty_\\%\n")
    for (chrom, ps_child), (h0_mismatches, h1_mismatches, n_variants, block_length) in n_sites.items():
        inherited_haplotype = "h0" if h0_mismatches < h1_mismatches else "h1" if h1_mismatches < h0_mismatches else ""
        certainty = 1 if ((h0_mismatches == 0 or h1_mismatches == 0) and not (h0_mismatches == 0 and h1_mismatches == 0)) else ""
        out.write(f"{chrom}\t{ps_child}\t{block_length}\t{n_variants}\t{h0_mismatches}\t{h1_mismatches}\t{inherited_haplotype}\t{certainty}\n")   
