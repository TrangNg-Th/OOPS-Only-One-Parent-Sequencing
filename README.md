# OOPS – Only One Parent Sequencing

A pipeline to find *de novo* mutation rates using only one parent's genetic
data.  It analyses sequencing reads (BAM) and variant calls (VCF) using
**samtools**, **bcftools**, and **whatshap**, orchestrated with **Snakemake**.

---

## Pipeline overview

| Step | Tool | Output |
|------|------|--------|
| Sort BAM | `samtools sort` | `results/sorted/*.sorted.bam` |
| Index BAM | `samtools index` | `results/sorted/*.sorted.bam.bai` |
| Alignment stats | `samtools flagstat` | `results/flagstat/*.flagstat.txt` |
| Variant calling | `bcftools mpileup \| bcftools call` | `results/vcf/*.raw.vcf.gz` |
| Variant filtering | `bcftools filter` | `results/vcf/*.filtered.vcf.gz` |
| Phase variants | `whatshap phase` | `results/phased/*.phased.vcf.gz` |
| Haplotag reads | `whatshap haplotag` | `results/haplotagged/*.haplotagged.bam` |
| Phasing stats | `whatshap stats` | `results/phase_stats/*.stats.tsv` |

---

## Requirements

[Conda](https://docs.conda.io/en/latest/) (or [Mamba](https://mamba.readthedocs.io/)) must be installed.

---

## Setup

```bash
# Create and activate the conda environment
conda env create -f environment.yml
conda activate oops-pipeline
```

---

## Usage

### 1. Prepare inputs

| Path | Description |
|------|-------------|
| `data/bam/{sample}.bam` | One unsorted or sorted BAM per sample |
| `data/ref/reference.fa` | Reference genome FASTA (indexed with `samtools faidx`) |

### 2. Configure the pipeline

Edit `config/config.yml`:

```yaml
samples:
  - sample1   # matches data/bam/sample1.bam
  - sample2

reference: "data/ref/reference.fa"
```

All other parameters (quality thresholds, ploidy, etc.) can be adjusted in the
same file.

### 3. Run the pipeline

```bash
# Dry-run (shows what will be executed without running anything)
snakemake --cores 8 --dry-run

# Full run
snakemake --cores 8
```

Use `--cores` to set the number of parallel threads/jobs.

---

## Output structure

```
results/
├── sorted/          # coordinate-sorted BAMs + indices
├── flagstat/        # samtools flagstat reports
├── vcf/             # raw and filtered VCFs
├── phased/          # whatshap-phased VCFs
├── haplotagged/     # read-haplotagged BAMs
├── phase_stats/     # whatshap phasing statistics (TSV)
└── logs/            # per-step log files
```
