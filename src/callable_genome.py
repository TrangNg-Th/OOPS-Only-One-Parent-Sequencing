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
          "<parent_id> <child_id> <GT_qual> <NV_quantile> "
          "<MM_diff_min> <WINDOW> <RECOUNT> [DNMC_FILE]")
    sys.exit(1)

in_file      = sys.argv[1]
out_dir      = sys.argv[2]
min_dp       = int(sys.argv[3])
max_dp       = int(sys.argv[4])
parent_id       = sys.argv[5]
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
print("Mother ID:", parent_id)
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
# Chromosome Exclusion (autosomes-only by default)
# =============================================================================
# Read from EXCLUDE_CHROMS env var (comma-separated). Defaults to "chrX,chrY,chrM".
# Matching is case-insensitive and tolerant of "chr"/"chrom" prefixes so users
# can pass chrX, chromX, X, etc
_excl_raw = os.environ.get("EXCLUDE_CHROMS", "chrX,chrY,chrM")

def _norm_chrom(name):
    s = str(name).strip().lower()
    if s.startswith("chrom"):
        s = s[5:]
    elif s.startswith("chr"):
        s = s[3:]
    return s

excluded = {_norm_chrom(x) for x in _excl_raw.split(",") if x.strip()}

if excluded:
    chrom_norm = df["chrom"].map(_norm_chrom)
    n_before = len(df)
    df = df[~chrom_norm.isin(excluded)]
    n_after = len(df)
    print(f"Excluding chromosomes: {sorted(excluded)} "
          f"(removed {n_before - n_after} SNPs, kept {n_after})")
else:
    print("No chromosome exclusion applied (EXCLUDE_CHROMS empty).")
    
    
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
summary_file = f"{out_dir}/{parent_id}_{child_id}_mismatch.{w}tsv"
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

# NOTE: use GroupBy.sample (not .apply(lambda x: x.sample(...))). Newer pandas
# excludes the grouping column ("PS_child") from the result of .apply(), which
# later crashes the consistency check with KeyError: 'PS_child'. GroupBy.sample
# keeps every column and is stable across pandas versions. Each PS block has >=1
# row by construction, so n=1 per group reproduces the old logic exactly.
sampled_df = (
    df.groupby("PS_child", group_keys=False, sort=False)
      .sample(n=1, random_state=RANDOM_SEED)
      .reset_index(drop=True)
)

print("Number of sampled denominator SNPs chosen from qualified blocks:", sampled_df.shape[0])

# total_sampled = sampled_df.shape[0]
# with open(callable_file, "a") as f:
#     f.write(f"total_sampled_snps\t{total_sampled}\n")
    
# =============================================================================
# Allele-balance filter on the SAMPLED set (AD[1:1]>5 && AD[1:0]>5, child)
# =============================================================================
# total_sampled (above) is the RANDOM denominator pool, taken before AD so the
# correction ratio's denominator stays the full sampled set. We now apply the
# child allele-balance filter to that sampled set, mirroring
# REfilter_dnm_candidates, so the BED written below contains only AD-passing
# sampled SNPs.
_AD_MIN = 5  # strictly greater-than, matching ">5" in the bcftools filter

for _c in ["ADref_child", "ADalt_child"]:
    if _c not in sampled_df.columns:
        print(f"ERROR: expected AD column '{_c}' not found in {in_file}. "
              "Re-extract _ps.tsv with the AD-aware query "
              "([%AD{0},][%AD{1},]) and update fix_PhaseSet.py headers.")
        sys.exit(1)

_adref_s = pd.to_numeric(sampled_df["ADref_child"], errors="coerce")
_adalt_s = pd.to_numeric(sampled_df["ADalt_child"], errors="coerce")

_n_before_ad = len(sampled_df)
sampled_df = sampled_df[(_adref_s > _AD_MIN) & (_adalt_s > _AD_MIN)].reset_index(drop=True)
print(f"Allele-balance filter on sampled set (child ADref>{_AD_MIN} & ADalt>{_AD_MIN}): "
      f"{_n_before_ad} -> {len(sampled_df)} SNPs")

total_sampled = sampled_df.shape[0]
with open(callable_file, "a") as f:
    f.write(f"total_sampled_snps\t{total_sampled}\n")
# =============================================================================
# Internal Consistency Check on Sampled SNPs
# =============================================================================
# Re-verify that every sampled SNP independently satisfies all filters this
# script applies, using the columns available in the input (_ps.tsv has no AD,
# so allele-balance is enforced later by REfilter_dnm_candidates, not here).
# This is a guard against format drift or filter-ordering regressions: if any
# sampled SNP fails, we fail loudly rather than write a tainted denominator.

def _check_sampled_consistency(s):
    qualified_blocks = set(summary_df["PS_child"].unique())
    problems = []

    for idx, row in s.iterrows():
        dp_mom   = int(row["DP_mom"])
        dp_child = int(row["DP_child"])
        gq_mom   = int(row["GQ_mom"])
        gq_child = int(row["GQ_child"])
        gt_mom   = row["GT_mom"]
        gt_child = row["GT_child"]
        chrom    = row["chrom"]
        cand     = row["cand_hapl"]
        ps       = row["PS_child"]

        # Depth bounds (both samples)
        if not (min_dp <= dp_mom <= max_dp):
            problems.append((idx, chrom, row["pos"], f"DP_mom={dp_mom} out of [{min_dp},{max_dp}]"))
        if not (min_dp <= dp_child <= max_dp):
            problems.append((idx, chrom, row["pos"], f"DP_child={dp_child} out of [{min_dp},{max_dp}]"))

        # Genotype quality (both samples)
        if gq_mom < GT_qual:
            problems.append((idx, chrom, row["pos"], f"GQ_mom={gq_mom} < {GT_qual}"))
        if gq_child < GT_qual:
            problems.append((idx, chrom, row["pos"], f"GQ_child={gq_child} < {GT_qual}"))

        # Biallelic 0/1-style GT (both samples)
        if not biallelic01.match(gt_mom):
            problems.append((idx, chrom, row["pos"], f"GT_mom={gt_mom} not biallelic 0/1-style"))
        if not biallelic01.match(gt_child):
            problems.append((idx, chrom, row["pos"], f"GT_child={gt_child} not biallelic 0/1-style"))

        # Autosome only
        if _norm_chrom(chrom) in excluded:
            problems.append((idx, chrom, row["pos"], f"chrom {chrom} is excluded"))

        # Membership in a qualified PS block
        if ps not in qualified_blocks:
            problems.append((idx, chrom, row["pos"], f"PS_child={ps} not in qualified blocks"))

        # cand_hapl / GT orientation (parent must be het; child phased to cand hap)
        parent_het = gt_mom in ("1/0", "0/1")
        if cand == "h0":
            ok = (gt_child == "1|0") and parent_het
        elif cand == "h1":
            ok = (gt_child == "0|1") and parent_het
        else:
            ok = False
        if not ok:
            problems.append((idx, chrom, row["pos"],
                             f"orientation fail: cand={cand} GT_child={gt_child} GT_mom={gt_mom}"))
            
        # Allele balance (child), mirrors REfilter_dnm_candidates AD[1:1]>5 & AD[1:0]>5
        adref_c = pd.to_numeric(pd.Series([row["ADref_child"]]), errors="coerce").iloc[0]
        adalt_c = pd.to_numeric(pd.Series([row["ADalt_child"]]), errors="coerce").iloc[0]
        if not (pd.notna(adref_c) and adref_c > _AD_MIN):
            problems.append((idx, chrom, row["pos"], f"ADref_child={row['ADref_child']} not > {_AD_MIN}"))
        if not (pd.notna(adalt_c) and adalt_c > _AD_MIN):
            problems.append((idx, chrom, row["pos"], f"ADalt_child={row['ADalt_child']} not > {_AD_MIN}"))
        

    return problems


_problems = _check_sampled_consistency(sampled_df)
if _problems:
    print("ERROR: sampled SNP consistency check FAILED for "
          f"{len(_problems)} violation(s):")
    for idx, chrom, pos, msg in _problems[:50]:
        print(f"  row {idx} {chrom}:{pos}  ->  {msg}")
    if len(_problems) > 50:
        print(f"  ... and {len(_problems) - 50} more")
    sys.exit(1)

print(f"Consistency check passed: all {sampled_df.shape[0]} sampled SNPs "
      "satisfy depth/GQ/GT/autosome/block/orientation/allele-balance filters.")

# =============================================================================
# Write BED Output
# =============================================================================
bed = pd.DataFrame({
    "chrom": sampled_df["chrom"],
    "start": sampled_df["pos"] - 1,
    "end": sampled_df["pos"]
}).drop_duplicates()

bed_file = f"{out_dir}/{parent_id}_{child_id}_hetc.{w}bed"
bed.to_csv(bed_file, sep="\t", index=False, header=False)

print("Qualified sampled snps written to bed file:", bed_file)

