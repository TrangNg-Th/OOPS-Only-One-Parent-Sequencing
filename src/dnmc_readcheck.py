import os
import sys
import pysam
from collections import Counter, defaultdict

os.chdir(os.getcwd())

# ============================================================
# Argument parsing
# ============================================================
if len(sys.argv) != 10:
    print(
        "Usage: python dnmc_readcheck.py "
        "<child_id> <child_bam> <input_bed> "
        "<min_base_qual> <min_mapping_qual> "
        "<window> <alt_read_count> <verbose> <label>"
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

if label == "dnmc":
    print("Processing DNM candidates...")
elif label == "hetc":
    print("Processing heterozygous candidates...")

# ============================================================
# Load candidate positions
# ============================================================
sites = []
with open(input_bed) as f:
    for line in f:
        chrom, start, end = line.rstrip().split("\t")[:3]
        sites.append((chrom, int(start), int(end)))

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

    # ========================================================
    # Store result
    # ========================================================
    ref_hp = next(iter(allele_hp[ref_allele]))
    alt_hp = next(iter(allele_hp[alt_allele]))
    
    
    final_results[(chrom, pos)] = {
        "REF": {
            "allele": ref_allele,
            "hp": next(iter(allele_hp[ref_allele])),
            "count": ref_count,
            "read_ids": allele_hp_reads[ref_allele][ref_hp],
        },
        "ALT": {
            "allele": alt_allele,
            "hp": next(iter(allele_hp[alt_allele])),
            "count": alt_count,
            "read_ids": allele_hp_reads[alt_allele][alt_hp],
        },
    }

bam.close()

# ============================================================
# Output summary
# ============================================================
print(f"Total candidates after read check: {len(final_results)}")
print("The candidates left are...")
for (chrom, pos), info in final_results.items():
    print(f"{chrom}:{pos} --> {info}")
    print()

# ============================================================
# Write BED (validated sites)
# ============================================================
out_bed = f"{child_id}_LR_validated_{label}.bed"

with open(out_bed, "w") as out:
    for (chrom, pos) in sorted(final_results):
        out.write(f"{chrom}\t{pos-1}\t{pos}\n")

print(f"Wrote {len(final_results)} sites to {out_bed}")

# ============================================================
# Write BED (windows for rephasing)
# ============================================================
size_kb = window // 1000
window_bed = f"{child_id}_{label}_plusminus{size_kb}kb.bed"

with open(window_bed, "w") as out:
    for (chrom, pos) in sorted(final_results):
        out.write(f"{chrom}\t{pos-window}\t{pos+window}\n")

print(f"Wrote {len(final_results)} sites to {window_bed}")


# ============================================================
# Write supporting read IDs for each haplotype / allele
# ============================================================
# Write supporting read IDs only when verbose == T and label == dnmc
if verbose and label == "dnmc":
    support_out = f"{child_id}_LR_validated_{label}_support_reads.tsv"

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