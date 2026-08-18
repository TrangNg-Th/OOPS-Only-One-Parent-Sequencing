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
# Optional keyword flag, pulled out of argv before the positional contract
# below is evaluated so the existing call sites keep working unchanged.
MAX_PARENT_AD = 0
if "--max-parent-ad" in sys.argv:
    _i = sys.argv.index("--max-parent-ad")
    MAX_PARENT_AD = int(sys.argv[_i + 1])
    del sys.argv[_i:_i + 2]

CLUSTER_WINDOW = 10000
if "--cluster-window" in sys.argv:
    _i = sys.argv.index("--cluster-window")
    CLUSTER_WINDOW = int(sys.argv[_i + 1])
    del sys.argv[_i:_i + 2]

# Reference-based low-complexity mask. A candidate sitting inside a homopolymer
# run is far more likely to be an alignment artifact than a real DNM: indels in
# such tracts are routinely represented as neighbouring substitutions.
MIN_HOMOPOLYMER = 6
if "--min-homopolymer" in sys.argv:
    _i = sys.argv.index("--min-homopolymer")
    MIN_HOMOPOLYMER = int(sys.argv[_i + 1])
    del sys.argv[_i:_i + 2]

REFERENCE = None
if "--reference" in sys.argv:
    _i = sys.argv.index("--reference")
    REFERENCE = sys.argv[_i + 1]
    del sys.argv[_i:_i + 2]

if len(sys.argv) not in (12,13,14): # one extra arg for passing the python script
    print("Usage: python count_mismatches.py "
          "<input_tsv> <out_dir> <min_dp> <max_dp> "
          "<parent_id> <child_id> <GT_qual> <NV_quantile> "
          "<MM_diff_min> <WINDOW> <RECOUNT> [DNMC_FILE]")

    print()
    print("Arguments:")
    print("  input_tsv     Input TSV file with genotype and phasing information")
    print("  out_dir       Output directory")
    print("  min_dp        Minimum read depth per site")
    print("  max_dp        Maximum read depth per site")
    print("  parent_id        Sample ID of the mother")
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
parent_id = sys.argv[5]
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
print("Mother ID:", parent_id)
print("Child ID:", child_id)
print("Genotype Quality threshold:", GT_qual)
print("Number of Variants Quantile threshold:", NV_QUANTILE)
print("Mismatch Difference Minimum:", MM_DIFF_MIN)
print("Max parent contradicting AD reads:", MAX_PARENT_AD)
print("Cluster-rejection window (bp):", CLUSTER_WINDOW)
print("Min homopolymer run to reject:", MIN_HOMOPOLYMER)
print("Reference for homopolymer mask:", REFERENCE)


print("--"*20)


# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
df = pd.read_csv(in_file, sep="\t", dtype=str)

# Add headers
# df.columns = ["chrom", "pos", "PS_mom", "PS_child", "GT_mom", "GT_child", "DP_mom", "DP_child", "GQ_mom", "GQ_child", "None"]


# ------------------------------------------------------------------
# Genotype sanity filtering
# ------------------------------------------------------------------
# Done before DP/GQ so that df_unfiltered below still has parseable genotypes.
biallelic01 = re.compile(r"^[01][\/|][01]$")
df = df[
    df["GT_mom"].str.match(biallelic01) &
    df["GT_child"].str.match(biallelic01)
]

df["pos"] = df["pos"].astype(int)

# ------------------------------------------------------------------
# Chromosome exclusion (must match the denominator in callable_genome.py)
# ------------------------------------------------------------------
# Part 4 excludes chrX/chrY/chrM from the callable genome. Without the same
# exclusion here the numerator could contain calls on chromosomes the
# denominator does not cover, which is a rate mismatch, not just noise.
_excl = {c.strip().lower() for c in
         os.environ.get("EXCLUDE_CHROMS", "chrX,chrY,chrM").split(",")
         if c.strip()}
if _excl:
    _n_before_excl = len(df)
    _keep = ~df["chrom"].str.lower().isin(_excl)
    df = df[_keep]
    print(f"Chromosome exclusion {sorted(_excl)}: "
          f"{_n_before_excl} -> {len(df)} sites")

# Snapshot before DP/GQ filtering. The cluster check further down needs to see
# the sites the depth/quality cuts are about to remove: a candidate whose
# neighbours were all filtered out looks like a lone mismatch when it is really
# one of a cluster, which is the signature of a local phasing artifact rather
# than a DNM.
df_unfiltered = df.copy()

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
# Parent AD veto (site-level)
# ------------------------------------------------------------------
# A homozygous parent call that still carries reads for the allele it does not
# claim is untrustworthy: a "mismatch" the child shows there may be an
# inherited variant the genotyper undercalled, not a DNM. Drop those sites
# only -- the PS block survives and simply loses a variant from n_variants.
# This is the VCF-AD counterpart to parent_readcheck.py, which reads the BAM
# with an MQ filter and so cannot see multi-mapping ALT support.
_ad_cols = ["ADref_mom", "ADalt_mom"]
if all(c in df.columns for c in _ad_cols):
    _adref_mom = pd.to_numeric(df["ADref_mom"], errors="coerce").fillna(0)
    _adalt_mom = pd.to_numeric(df["ADalt_mom"], errors="coerce").fillna(0)
    _gt_mom = df["GT_mom"].str.replace("|", "/", regex=False)

    parent_ad_bad = (
        ((_gt_mom == "0/0") & (_adalt_mom > MAX_PARENT_AD)) |
        ((_gt_mom == "1/1") & (_adref_mom > MAX_PARENT_AD))
    )
    print(f"Parent AD veto (hom parent with > {MAX_PARENT_AD} contradicting "
          f"reads): {len(df)} -> {len(df) - int(parent_ad_bad.sum())} sites")
    df = df[~parent_ad_bad]
else:
    print(f"WARNING: {_ad_cols} not found in {in_file}; skipping parent AD veto. "
          "Re-extract _ps.tsv with the AD-aware query ([%AD{0},][%AD{1},]).")

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
quantiles = n_variants_series.quantile([0, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99])

print("n_variants quantiles:")
print(quantiles)
print("Writing file *mismatch.tsv")
summary_file = f"{out_dir}/{parent_id}_{child_id}_mismatch.{w}tsv"
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
else:
    summary_df = summary_df[summary_df["mismatch_difference"] < 0.5] 
    # when recounting, we allow for more similar mismatch fractions between the two blocks, as the window is smaller and thus the number of variants is smaller, which can lead to more stochasticity in the mismatch fractions

# ------------------------------------------------------------------
# Threshold 4: Only for recounting DNMs, we require that the max_mm has to be at least 50% 
# of the number of variants in the block, to ensure that there is a strong signal of
# mismatch in at least one block (since we use smaller blocks)

# if RECOUNT == "T":
#     summary_df = summary_df[
#         (summary_df[["n_mismatches_h0", "n_mismatches_h1"]].max(axis=1) >= 0.20 * summary_df["n_variants"])
#     ]
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

# ------------------------------------------------------------------
# Cluster index (built on the pre-DP/GQ snapshot)
# ------------------------------------------------------------------
# For every PS block, record the positions where each child haplotype
# mismatches the parent, using sites that have NOT yet been through the
# depth/quality cuts. A genuine DNM is a solitary event; if the same
# haplotype mismatches again within CLUSTER_WINDOW bp, the block has a local
# phasing/mapping problem and the candidate is not trustworthy.
raw_mm_by_block = {}
if CLUSTER_WINDOW > 0:
    for (chrom, ps), g in df_unfiltered.groupby(["chrom", "PS_child"]):
        v = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
        h0, h1 = v[0].values, v[1].values
        mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

        mask_h0 = (mom == h0[:, None]).any(axis=1)
        mask_h1 = (mom == h1[:, None]).any(axis=1)

        pos_arr = g["pos"].values
        raw_mm_by_block[(chrom, ps)] = (pos_arr[~mask_h0], pos_arr[~mask_h1])


def has_cluster_neighbour(chrom, ps, pos, hap_idx):
    """True if the same haplotype mismatches again within CLUSTER_WINDOW bp.

    hap_idx is 0 or 1 and indexes the child haplotype carrying the candidate.
    Phasing is consistent within a PS block, so the index is comparable across
    all sites of that block.
    """
    if CLUSTER_WINDOW <= 0:
        return False
    entry = raw_mm_by_block.get((chrom, ps))
    if entry is None:
        return False
    others = entry[hap_idx]
    near = others[(np.abs(others - pos) <= CLUSTER_WINDOW) & (others != pos)]
    return near.size > 0


_fasta = None
if MIN_HOMOPOLYMER > 0 and REFERENCE:
    try:
        import pysam
        _fasta = pysam.FastaFile(REFERENCE)
    except Exception as e:
        print(f"WARNING: cannot open reference {REFERENCE} ({e}); "
              "homopolymer mask disabled")
        _fasta = None
elif MIN_HOMOPOLYMER > 0:
    print("WARNING: --min-homopolymer set but no --reference; mask disabled")


def in_homopolymer(chrom, pos, flank=60):
    """True if pos sits in (or immediately beside) a long homopolymer run.

    A de novo SNV is a single base change; a run of identical reference bases
    around it means the aligner had a length-ambiguous tract to place reads in,
    which is where spurious neighbouring substitutions come from.
    """
    if _fasta is None or MIN_HOMOPOLYMER <= 0:
        return False
    start = max(0, pos - 1 - flank)
    try:
        seq = _fasta.fetch(chrom, start, pos + flank).upper()
    except Exception:
        return False
    if not seq:
        return False
    idx = (pos - 1) - start          # candidate offset within seq
    # Longest run overlapping the candidate base or either neighbour.
    for probe in (idx - 1, idx, idx + 1):
        if probe < 0 or probe >= len(seq):
            continue
        base = seq[probe]
        if base not in "ACGT":
            continue
        lo = hi = probe
        while lo > 0 and seq[lo - 1] == base:
            lo -= 1
        while hi < len(seq) - 1 and seq[hi + 1] == base:
            hi += 1
        if (hi - lo + 1) >= MIN_HOMOPOLYMER:
            return True
    return False


# Save results
dnm_records = []
n_cluster_rejected = 0
n_homopolymer_rejected = 0

if RECOUNT != "T":
    for (chrom, ps), g in merged.groupby(["chrom", "PS_child"]):
        v = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
        h0, h1 = v[0].values, v[1].values
        mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

        mask_h0 = (mom == h0[:, None]).any(axis=1)
        mask_h1 = (mom == h1[:, None]).any(axis=1)

        n0, n1 = np.sum(~mask_h0), np.sum(~mask_h1)
        if n0 == 1 and n1 > 1:
            pos = g.loc[~mask_h0, "pos"].iloc[0]
            hap_idx = 0
        elif n1 == 1 and n0 > 1:
            pos = g.loc[~mask_h1, "pos"].iloc[0]
            hap_idx = 1
        else:
            continue

        if has_cluster_neighbour(chrom, ps, pos, hap_idx):
            n_cluster_rejected += 1
            continue

        if in_homopolymer(chrom, pos):
            n_homopolymer_rejected += 1
            continue

        dnm_records.append({
            "chrom": chrom,
            "pos": pos,
            "PS_child": ps
        })
elif RECOUNT == "T":
    print(f"Recounting candidate DNMs from file {DNMC_FILE} with window {WINDOW}bp")
    
    for (chrom, ps), g in merged.groupby(["chrom", "PS_child"]):
        
        v = g["GT_child"].str.split(r"[\/|]", expand=True).astype(int)
        h0, h1 = v[0].values, v[1].values
        mom = g["GT_mom"].str.split(r"[\/|]", expand=True).astype(int).values

        mask_h0 = (mom == h0[:, None]).any(axis=1)
        mask_h1 = (mom == h1[:, None]).any(axis=1)

        # The logic here is as follows:
        # - If one block has exactly one mismatch and the other block has more than one mismatch
        if (np.sum(~mask_h0) == 1) and (np.sum(~mask_h1) > 1):
            pos = g.loc[~mask_h0, "pos"].iloc[0]
            hap_idx = 0

        # - If one block has exactly one mismatch and the other block has more than one mismatch
        elif (np.sum(~mask_h1) == 1) and (np.sum(~mask_h0) > 1):
            pos = g.loc[~mask_h1, "pos"].iloc[0]
            hap_idx = 1

        else:
            continue


        if (chrom, pos) not in dnmc_set:
            continue

        if has_cluster_neighbour(chrom, ps, pos, hap_idx):
            n_cluster_rejected += 1
            continue

       
        if in_homopolymer(chrom, pos):
            n_homopolymer_rejected += 1
            continue

        dnm_records.append({
            "chrom": chrom,
            "pos": pos,
            "PS_child": ps
        })
print(f"Homopolymer rejection (run >= {MIN_HOMOPOLYMER}bp at candidate): "
      f"{n_homopolymer_rejected} candidates dropped")
print(f"Cluster rejection (same haplotype mismatching again within "
      f"{CLUSTER_WINDOW}bp, pre-DP/GQ): {n_cluster_rejected} candidates dropped")

if not dnm_records:
    print("No candidate DNMs found after recounting with the specified criteria.")
else:
    dnm_df = pd.DataFrame(dnm_records)
    dnm_df.to_csv(f"{out_dir}/{parent_id}_{child_id}_dnmc.{w}tsv", sep="\t", index=False)
    print("Writing file *dnmc.tsv")
    # BED
    bed = pd.DataFrame({
        "chrom": dnm_df["chrom"],
        "start": dnm_df["pos"]-1,
        "end": dnm_df["pos"]
    }).drop_duplicates()

    print("Writing file *dnmc.bed")
    bed.to_csv(f"{out_dir}/{parent_id}_{child_id}_dnmc.{w}bed", sep="\t", index=False, header=False)


