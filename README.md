
# One-Parent Sequencing Mutation Rate Pipeline

> Estimating de novo mutation rates using child long reads, child short reads, and a single parent short-read dataset.

---

## Overview

This pipeline estimates the **de novo mutation (DNM) rate** in a child using:

* Child **long reads** (HiFi BAM)
* Child **short reads** (Illumina VCF )
* One **parent short-read** dataset (Illumina VCF)

Unlike classical trio-based pipelines, this workflow is designed for scenarios where **only one parent is sequenced**, leveraging:

* Long-read–guided phasing
* Haplotype block analysis
* Local rephasing around candidate DNMs

Final mutation rate:

```
Mutation Rate = (# Validated DNMs) / (Callable genome size)
```

---

# Pipeline Architecture

```
PART 1  →  Data preparation
PART 2  →  Child phasing using long reads (HiFi Bam)
PART 2b →  Initial DNM detection
PART 3  →  Local rephasing around DNMs
PART 3b →  Refined DNM detection 
PART 4  →  Callable genome estimation & mutation rate
```

Each stage can be run independently.

---

# Repository Structure

```
.
├── pipeline.sh
├── README.md
└── variants/
    └── src/
        ├── count_mismatches.py
        ├── callable_genome.py
        ├── dnmc_readcheck.py
        └── fix_PhaseSet.py
```

---

# Project Directory Layout (Runtime)

The runtime directory is defined via:

```
--prj-dir
```

Example:

```
/project/mutation_rate/
```

*** Expected structure:***

```
PRJ_DIR/
│
├── reference/
│   └── chm13v2.0_maskedY_rCRS.fa
│
├── hifi/
│   └── SAMPLE_CHILD.CHM13.haplotagged.bam
│
└── variants/
    └── small_variants/
        ├── illumina-dragen/
        ├── hifi/
        └── SAMPLE_CHILD_phasedvcf/
```

---

# Directory Details

## 1. `reference/`

Contains reference genome used for:

* Phasing
* VCF normalization
* Long-read validation

```
chm13v2.0_maskedY_rCRS.fa
```

---

## 2. `hifi/`

Child long-read BAM:

```
SAMPLE_CHILD.CHM13.haplotagged.bam
SAMPLE_CHILD.CHM13.haplotagged.bam.bai
```

Used for:

* Whatshap phasing
* Long-read validation of DNM candidates

---

## 3. `variants/small_variants/illumina-dragen/`

Contains short-read VCFs.

### After splitting:

```
SAMPLE_PARENT.CHM13.illumina-dragen.oa.vcf.gz
SAMPLE_CHILD.CHM13.illumina-dragen.oa.vcf.gz
```

### After removing phasing artifacts:

```
*.unphased.noPS.vcf.gz
```

---

## 4. `variants/small_variants/SAMPLE_CHILD_phasedvcf/`

Main analysis directory:

```
SAMPLE_CHILD_phasedvcf/
│
├── SAMPLE_CHILD.illumVCF_hifiOnly.phased.vcf.gz
├── SAMPLE_CHILD.illumVCF_hifiOnly.stats.tsv
├── SAMPLE_CHILD.illumVCF_hifiOnly.blocks.tsv
│
├── merged/
├── HTblocks/
├── mismatch_analysis/
```

---

### `merged/`

Merged parent + child VCF:

```
SAMPLE_PARENT_SAMPLE_CHILD.merged.vcf.gz
```

Used for haplotype comparison.

---

### `HTblocks/`

Extracted haplotype blocks:

```
*_ps.tsv
```

Contains:

| Column | Description      |
| ------ | ---------------- |
| CHROM  | Chromosome       |
| POS    | Position         |
| PS     | Phase set ID     |
| GT     | Genotype         |
| DP     | Depth            |
| GQ     | Genotype quality |

---

### `mismatch_analysis/`

Core DNM detection outputs:

```
*_dnmc.tsv
*_dnmc.bed
*_dnmc.vcf.gz
*_LR_validated.tsv
```

#### `_dnmc.tsv`

Candidate DNM list.

#### `_dnmc.bed`

BED file of DNM candidate positions.

#### `_dnmc.vcf.gz`

Subset VCF restricted to DNM candidates.

#### `_LR_validated.tsv`

Long-read validation results.

---

### `denum_calcul/`

Callable genome analysis:

```
*_hetc.tsv
*_hetc.bed
*_hetc.vcf.gz
```

Used to compute the denominator of mutation rate.

---

# How the Method Works

### Step 1 — Phase the child

Child short reads are phased using long reads via **Whatshap**.

### Step 2 — Merge with one parent

Phased child VCF is merged with unphased parent VCF.

### Step 3 — Haplotype mismatch detection

Within each phase block:

* Parent = homozygous reference
* Child = heterozygous on a specific haplotype

These are candidate DNMs.

### Step 4 — Long-read validation

Candidates must:

* Have ALT support in long reads
* Pass base quality threshold
* Pass mapping quality threshold

### Step 5 — Local rephasing

Rephase ±window around each candidate to confirm haplotype consistency.

### Step 6 — Callable genome estimation

Apply same filters to randomly sampled heterozygous SNPs to estimate:

```
Callable genome size
```

---

# Running the Pipeline

## Required arguments

```
--prj-dir
--sample-child
--sample-parent
```

## Examples of command
Full commands to run all parts of the pipeline
```bash
./main.sh \
  --prj-dir /project/mutation_rate \
  --sample-child NA12879 \
  --sample-parent NA12878 \
  --all
```

```bash 
./main.sh \
  --all \
  --prj-dir /N/project/mutation_rate_Mmulatta/platinum-ped-data/aws-data \
  --sample-child NA12879 \
  --sample-parent NA12878 \
  --cpus 8 \
  --time 4:00:00 \
  --min-rdepth 15 \
  --max-rdepth 50 \
  --gt-qual 30 \
  --nv-quantile 0.75 \
  --mm-diff-min 0.1 \
  --min-base-qual 20 \
  --min-map-qual 20 \
  --window 20000 \
  --alt-read-count 8 \
  --verbose T
```

## Run specific parts

Run only data preparation + phasing
```bash 
./main.sh \
  --part 1 2 \
  --prj-dir /N/project/mutation_rate_Mmulatta/platinum-ped-data/aws-data \
  --sample-child NA12879 \
  --sample-parent NA12878
```

Run only the DNM detection 
```bash 
./pipeline.sh \
  --part 2b \
  --prj-dir /N/project/mutation_rate_Mmulatta/platinum-ped-data/aws-data \
  --sample-child NA12879 \
  --sample-parent NA12878
```

---

# Software Requirements

* SLURM
* bcftools
* whatshap
* Python 3
* AWS CLI (optional)
* Conda

---

# Outputs

After full run:

* Final validated DNM list
* Callable genome size
* Mutation rate estimate

---

# Scientific Assumptions

* Accurate long-read phasing
* Parent truly homozygous reference at DNM sites
* Phase blocks reflect true haplotypes
* Filters applied identically to numerator and denominator

---

# Reproducibility

* Modular execution (`--part`)
* Deterministic filtering
* Explicit denominator calculation
* Long-read confirmation stage

---

# Citation

If using this pipeline in a publication, please cite:

> [Nguyen], One-Parent Sequencing Mutation Rate Pipeline, GitHub repository.


---

# License
None

---

# Future Improvements

* Multi-parent extension
* Automatic Slurm dependency chaining
* Containerization (Docker/Singularity)
* Workflow manager integration (Nextflow/Snakemake)

---
