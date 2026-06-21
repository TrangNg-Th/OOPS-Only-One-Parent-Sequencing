#!/usr/bin/env python3

import os
import sys
import pysam
from collections import Counter, defaultdict

# ============================================================
# Argument parsing
# ============================================================
# Positional args (1-indexed as seen on the command line):
#   1  child_id
#   2  child_bam
#   3  input_bed
#   4  min_base_qual
#   5  min_mapping_qual
#   6  window
#   7  alt_read_count
#   8  verbose            (T/F)
#   9  label              (e.g. dnmc, hetc, dnmc_rephase, hetc_rephase)
#   10 total_read_min     (minimum ref+alt support; optional, may be absent)
#   11 out_dir            (directory to write outputs into; optional)
#
# out_dir is OPTIONAL and defaults to the current working directory, so older
# callers that pass only 10 args keep working unchanged. New callers pass an
# explicit directory so outputs land next to the rest of that step's results
# instead of in whatever directory the script happened to be launched from.
if len(sys.argv) < 11 or len(sys.argv) > 12:
    print(
        "Usage: python dnmc_readcheck.py "
        "<child_id> <child_bam> <input_bed> "
        "<min_base_qual> <min_mapping_qual> "
        "<window> <alt_read_count> <verbose> <label> <total_read_min> "
        "[out_dir]"
    )
    sys.exit(1)

child_id = sys.argv[1]
child_bam = sys.argv[2]
input_bed = sys.argv[3]
min_base_qual = int(sys.argv[4])
min_mapping_qual = int(sys.argv[5])
window = int(sys.argv[6])
alt_read_count = int(sys.argv[7])
verbose = sys.argv[8].upper() == "T"
label = sys.argv[9].lower()
total_read_min = sys.argv[10] if len(sys.argv) > 10 else None

# Output directory: arg 11 if given, else CWD (back-compatible default).
out_dir = sys.argv[11] if len(sys.argv) > 11 else os.getcwd()
os.makedirs(out_dir, exist_ok=True)

# ------------------------------------------------------------
# Label handling
# ------------------------------------------------------------
# We accept the original labels ("dnmc", "hetc") AND their rephase
# variants ("dnmc_rephase", "hetc_rephase"). The "base kind" drives
# behaviour (banner text + whether verbose support reads are written),
# while the full label is preserved verbatim in output filenames so
# that set-A and set-B (rephase) runs never collide.
base_kind = label
for suffix in ("_rephase",):
    if base_kind.endswith(suffix):
        base_kind = base_kind[: -len(suffix)]
        break

if base_kind == "dnmc":
    print(f"Processing DNM candidates... (label='{label}')")
elif base_kind == "hetc":
    print(f"Processing heterozygous candidates... (label='{label}')")
else:
    print(f"Processing candidates with label='{label}'")

print(f"Output directory: {out_dir}")

# ============================================================
# Load candidate positions
# ============================================================
sites = []
if not os.path.exists(input_bed):
    print(f"ERROR: input BED not found: {input_bed}")
    sys.exit(1)

with open(input_bed) as f:
    for line in f:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        # tolerate a header line if one ever sneaks in
        if line.startswith("#") or line.lower().startswith("chrom\t"):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        chrom, start, end = fields[:3]
        sites.append((chrom, int(start), int(end)))

if len(sites) == 0:
    print(f"WARNING: no candidate sites found in {input_bed}")

bam = pysam.AlignmentFile(child_bam, "rb")

final_results = {}

# ============================================================
# Main loop: evaluate each site
# ============================================================
for chrom, start, end in sites:
    pos = end  # 1-based coordinate

    # allele -> Counter(haplotype -> read count)
    allele_hp = defaultdict(Counter)
    allele_hp_reads = defaultdict(lambda: defaultdict(list))  # allele -> hp -> list of read names

    for col in bam.pileup(
        chrom,
        start,
        end,
        min_base_quality=min_base_qual,
        min_mapping_quality=min_mapping_qual,
        truncate=True,
    ):
        if col.pos != pos - 1:
            continue

        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue

            read = pr.alignment
            base = read.query_sequence[pr.query_position].upper()

            if base not in "ACGT":
                continue

            hp = read.get_tag("HP") if read.has_tag("HP") else None
            allele_hp[base][hp] += 1
            allele_hp_reads[base][hp].append(read.query_name)

    ## Debugging output
    if verbose:
        print(f"Processing {chrom}:{pos}")
        for allele, hp_counts in allele_hp.items():
            print(f"{allele} -> {hp_counts}")
        print("=" * 40)

    # ========================================================
    # Filtering logic (biological constraints)
    # ========================================================

    # Must have exactly 2 alleles
    if len(allele_hp) != 2:
        continue

    # Each allele must map to exactly one haplotype
    if not all(len(hp_counts) == 1 for hp_counts in allele_hp.values()):
        continue

    # Alleles must be on different haplotypes
    haplotypes = [next(iter(hp_counts)) for hp_counts in allele_hp.values()]
    if len(set(haplotypes)) != 2:
        continue

    # ========================================================
    # Compute allele counts
    # ========================================================
    allele_totals = {
        allele: sum(hp_counts.values())
        for allele, hp_counts in allele_hp.items()
    }

    alleles = list(allele_totals.keys())
    if len(alleles) != 2:
        continue

    a1, a2 = alleles
    c1, c2 = allele_totals[a1], allele_totals[a2]

    # Reject ambiguous tie cases
    if c1 == c2:
        continue

    # Assign REF = more supported allele, ALT = less supported
    if c1 > c2:
        ref_allele, alt_allele = a1, a2
        ref_count, alt_count = c1, c2
    else:
        ref_allele, alt_allele = a2, a1
        ref_count, alt_count = c2, c1

    # Require minimum support for ALT allele
    if alt_count < alt_read_count:
        continue

    # Require minimum total read support (if specified)
    if total_read_min is not None:
        total_count = ref_count + alt_count
        if total_count < int(total_read_min) - 1:
            continue

    # ========================================================
    # Store result
    # ========================================================
    ref_hp = next(iter(allele_hp[ref_allele]))
    alt_hp = next(iter(allele_hp[alt_allele]))

    final_results[(chrom, pos)] = {
        "REF": {
            "allele": ref_allele,
            "hp": ref_hp,
            "count": ref_count,
            "read_ids": allele_hp_reads[ref_allele][ref_hp],
        },
        "ALT": {
            "allele": alt_allele,
            "hp": alt_hp,
            "count": alt_count,
            "read_ids": allele_hp_reads[alt_allele][alt_hp],
        },
    }

bam.close()

# ============================================================
# Output summary
# ============================================================
print(f"Total candidates after read check: {len(final_results)}")
if verbose:
    print("The candidates left are...")
    for (chrom, pos), info in final_results.items():
        print(f"{chrom}:{pos} --> {info}")
        print()

# ============================================================
# Write BED (validated sites)
# ============================================================
out_bed = os.path.join(out_dir, f"{child_id}_LR_validated_{label}.bed")

with open(out_bed, "w") as out:
    for (chrom, pos) in sorted(final_results):
        out.write(f"{chrom}\t{pos-1}\t{pos}\n")

print(f"Wrote {len(final_results)} sites to {out_bed}")

# ============================================================
# Write BED (windows for rephasing)
# ============================================================
size_kb = window // 1000
window_bed = os.path.join(out_dir, f"{child_id}_{label}_plusminus{size_kb}kb.bed")

with open(window_bed, "w") as out:
    for (chrom, pos) in sorted(final_results):
        out.write(f"{chrom}\t{pos-window}\t{pos+window}\n")

print(f"Wrote {len(final_results)} sites to {window_bed}")


# ============================================================
# Write supporting read IDs for each haplotype / allele
# ============================================================
# Write supporting read IDs only when verbose == T and this is a DNM run
# (base_kind == "dnmc"), so both the original "dnmc" and the rephase
# "dnmc_rephase" runs produce the support file.
if verbose and base_kind == "dnmc":
    support_out = os.path.join(out_dir, f"{child_id}_LR_validated_{label}_support_reads.tsv")

    with open(support_out, "w") as out:
        out.write("chrom\tpos\ttype\tallele\thaplotype\tread_count\tread_ids\n")

        for (chrom, pos), info in sorted(final_results.items()):
            for allele_type in ["REF", "ALT"]:
                allele = info[allele_type]["allele"]
                hp = info[allele_type]["hp"]
                count = info[allele_type]["count"]
                read_ids = info[allele_type]["read_ids"]

                out.write(
                    f"{chrom}\t{pos}\t{allele_type}\t{allele}\t{hp}\t{count}\t{','.join(read_ids)}\n"
                )

    print(f"Wrote supporting read IDs to {support_out}")