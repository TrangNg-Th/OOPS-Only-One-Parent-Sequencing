#!/usr/bin/env python3
"""
Direct parent BAM validation for DNM candidates.

For each candidate site:
  - look up REF/ALT from the merged VCF
  - pileup the parent BAM and count reads supporting each allele
  - reject the candidate if the parent has > max_alt_reads ALT-supporting reads

Catches false-negative parental 0/0 calls (where DRAGEN dropped ALT reads
due to nearby variants, soft clipping, or local realignment issues), which
the bcftools VCF-level filter structurally cannot see.
"""
import os
import sys
import pysam
from collections import Counter

if len(sys.argv) != 10:
    print(
        "Usage: python parent_readcheck.py "
        "<parent_id> <parent_bam> <merged_vcf> "
        "<input_bed> <output_bed> "
        "<min_base_qual> <min_mapping_qual> <max_alt_reads> <out_dir>"
    )
    sys.exit(1)

parent_id         = sys.argv[1]
parent_bam_path   = sys.argv[2]
merged_vcf_path   = sys.argv[3]
input_bed         = sys.argv[4]
output_bed        = sys.argv[5]
min_base_qual     = int(sys.argv[6])
min_mapping_qual  = int(sys.argv[7])
max_alt_reads     = int(sys.argv[8])
out_dir           = sys.argv[9]

os.makedirs(out_dir, exist_ok=True)

print("-" * 60)
print(f"[parent-check] Parent ID          : {parent_id}")
print(f"[parent-check] Parent BAM         : {parent_bam_path}")
print(f"[parent-check] Merged VCF         : {merged_vcf_path}")
print(f"[parent-check] Input BED          : {input_bed}")
print(f"[parent-check] Output BED         : {output_bed}")
print(f"[parent-check] Min base qual      : {min_base_qual}")
print(f"[parent-check] Min mapping qual   : {min_mapping_qual}")
print(f"[parent-check] Max parental ALT   : {max_alt_reads}")
print("-" * 60)

# Load candidate sites
sites = []
with open(input_bed) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3 or parts[0].startswith("#"):
            continue
        sites.append((parts[0], int(parts[1]), int(parts[2])))

if not sites:
    print("[parent-check] No candidate sites in input BED. Writing empty output.")
    open(output_bed, "w").close()
    sys.exit(0)

vcf = pysam.VariantFile(merged_vcf_path)
bam = pysam.AlignmentFile(parent_bam_path, "rb")

passed   = []
rejected = []

for chrom, start, end in sites:
    pos = end  # 1-based

    # --- look up REF/ALT in merged VCF ----------------------------------
    ref, alt = None, None
    for rec in vcf.fetch(chrom, pos - 1, pos):
        if rec.pos != pos:
            continue
        if rec.ref is None or rec.alts is None:
            continue
        # Restrict to biallelic SNVs (matches upstream `bcftools view -m2 -M2 -v snps`)
        if len(rec.ref) == 1 and len(rec.alts) == 1 and len(rec.alts[0]) == 1:
            ref, alt = rec.ref.upper(), rec.alts[0].upper()
        break

    if ref is None or alt is None:
        rejected.append((chrom, start, end, "NA", "NA", 0, 0, "no_biallelic_snv_in_vcf"))
        continue

    # --- pileup parent BAM ----------------------------------------------
    counts = Counter()
    for col in bam.pileup(
        chrom, pos - 1, pos,
        min_base_quality=min_base_qual,
        min_mapping_quality=min_mapping_qual,
        truncate=True,
    ):
        if col.pos != pos - 1:
            continue
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue
            base = pr.alignment.query_sequence[pr.query_position].upper()
            if base in "ACGT":
                counts[base] += 1

    ref_n = counts.get(ref, 0)
    alt_n = counts.get(alt, 0)
    other = sum(v for b, v in counts.items() if b not in (ref, alt))

    print(f"[parent-check] {chrom}:{pos}  REF={ref}({ref_n})  ALT={alt}({alt_n})  other={other}", end="")

    if alt_n > max_alt_reads:
        print("  -> REJECT (parent carries ALT)")
        rejected.append((chrom, start, end, ref, alt, ref_n, alt_n, f"parent_alt={alt_n}"))
    else:
        print("  -> PASS")
        passed.append((chrom, start, end, ref, alt, ref_n, alt_n))

vcf.close()
bam.close()

# --- Write filtered BED -------------------------------------------------
with open(output_bed, "w") as out:
    for chrom, start, end, *_ in passed:
        out.write(f"{chrom}\t{start}\t{end}\n")

# --- Audit TSV (everything we inspected, PASS + REJECT, with counts) ----
audit_tsv = os.path.join(out_dir, f"{parent_id}_parent_BAM_check.tsv")
with open(audit_tsv, "w") as out:
    out.write("chrom\tpos\tref\talt\tparent_ref_reads\tparent_alt_reads\tdecision\treason\n")
    for chrom, start, end, ref, alt, ref_n, alt_n in passed:
        out.write(f"{chrom}\t{end}\t{ref}\t{alt}\t{ref_n}\t{alt_n}\tPASS\t-\n")
    for chrom, start, end, ref, alt, ref_n, alt_n, reason in rejected:
        out.write(f"{chrom}\t{end}\t{ref}\t{alt}\t{ref_n}\t{alt_n}\tREJECT\t{reason}\n")

print("-" * 60)
print(f"[parent-check] Sites checked : {len(sites)}")
print(f"[parent-check] Passed        : {len(passed)}")
print(f"[parent-check] Rejected      : {len(rejected)}")
print(f"[parent-check] Output BED    : {output_bed}")
print(f"[parent-check] Audit TSV     : {audit_tsv}")
print("-" * 60)