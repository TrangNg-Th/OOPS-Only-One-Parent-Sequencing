import os
import sys
import pysam  ## included in whatshap-env
from collections import Counter, defaultdict


# set working directory


## ====================================================
## older scripts
## =======================================================
os.chdir(os.getcwd())
# child_bam = "../../hifi/NA12879.CHM13.haplotagged.bam"
# input_dnmc_bed = "../small_variants/NA12879_phasedvcf/mismatch_analysis/NA12878_NA12879_dnmc.bed"

# min_base_qual = 20
# min_mapping_qual = 20
## =====================================================
if len(sys.argv) != 9:
    print(
        "Usage: python dnmc_readcheck.py ",
        "<child_bam> <input_dnmc_bed> <min_base_qual> <min_mapping_qual> "
    )
    
    print()
    print("Arguments:")
    print("  child_id           Input ID of the child")
    print("  child_bam          Input long reads bam file for the child")
    print("  input_dnmc_bed     Input DNM candidate positions")
    print("  min_base_qual      Input minimum base quality threshold (ex. 20)")
    print("  min_mapping_qual   Input minimum mapping quality threshold (ex. 20)")
    print("  window             Window to create BED file of pos +/- window (ex. 10000)")
    print("  alt_read_count     Number of reads of a HP supporting the alternative allele")
    print("  verbose (T or F)   If T, print out all the alleles counted")
    
    


## =============================================================
## Parse arguments
## =============================================================
child_id = sys.argv[1]
child_bam = sys.argv[2]
input_dnmc_bed = sys.argv[3]
min_base_qual = int(sys.argv[4])
min_mapping_qual = int(sys.argv[5])
window = int(sys.argv[6])
alt_read_count = int(sys.argv[7])
verbose=sys.argv[8]

if verbose == "T":
    verbose = True
else:
    verbose = False

## ===========================================================


# Read DNM candidate BED
dnmc_sites = []
with open(input_dnmc_bed) as f:
    for line in f:
        chrom, start, end = line.rstrip().split("\t")[:3]
        dnmc_sites.append((chrom, int(start), int(end)))

bam = pysam.AlignmentFile(child_bam, "rb")

final_dnmc_result = {}

for chrom, start, end in dnmc_sites:
    pos = end  # 1-based
    

    # allele -> Counter(HP -> read_count)
    allele_hp = defaultdict(Counter)

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
            
    if verbose: 
        print(f"Processing chrom {chrom} @ {end}")
        for k, v in allele_hp.items():
            print(f"{k} --> {v}")
            print()
        print("=====================================")
    
    # ---- strict phasing criterion ----
    if len(allele_hp) != 2:
        continue

    # each allele must be on exactly one haplotype
    if not all(len(hp_counts) == 1 for hp_counts in allele_hp.values()):
        continue

    # alleles must be on different haplotypes
    haplotypes = [next(iter(hp_counts)) for hp_counts in allele_hp.values()]
    if len(set(haplotypes)) != 2:
        continue

    # KEEP: store full read-count info
     # ---- determine REF / ALT by read count ----
    allele_totals = {
        allele: sum(hp_counts.values())
        for allele, hp_counts in allele_hp.items()
    }

    # ALT = allele with fewer reads
    alt_allele = min(allele_totals, key=allele_totals.get)
    ref_allele = max(allele_totals, key=allele_totals.get)

    alt_count = allele_totals[alt_allele]

    # ---- ALT read count filter ----
    if alt_count <= alt_read_count:  # at least n read per alternative allele in the long read
        continue

    # KEEP: store tagged result
    final_dnmc_result[(chrom, pos)] = {
        "REF": {
            "allele": ref_allele,
            "hp": next(iter(allele_hp[ref_allele])),
            "count": allele_totals[ref_allele],
        },
        "ALT": {
            "allele": alt_allele,
            "hp": next(iter(allele_hp[alt_allele])),
            "count": alt_count,
        },
    }

bam.close()

print("The candidates left are...")
for k, v in final_dnmc_result.items():
    print(f"{k[0]}:{k[1]} --> {v}")
    print()
    
# ## Create bed file for the candidates
# ## ===========================================================
out_bed = f"{child_id}_LR_validated_dnmc.bed"

with open(out_bed, "w") as out:
    for (chrom, pos) in sorted(final_dnmc_result):
        out.write(f"{chrom}\t{pos-1}\t{pos}\n")
print(f"Wrote {len(final_dnmc_result)} sites to {out_bed}")

## Create bed file for rephasing
# ## ===========================================================
size=window/1000
out_bed = f"{child_id}_dnm_plusminus{int(size)}kb.bed"

with open(out_bed, "w") as out:
    for (chrom, pos) in sorted(final_dnmc_result):
        out.write(f"{chrom}\t{pos-window}\t{pos+window}\n")
print(f"Wrote {len(final_dnmc_result)} sites to {out_bed}")


# Write out the bed file of the candidate
# print(final_dnmc_result)