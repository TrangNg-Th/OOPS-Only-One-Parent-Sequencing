Here’s a **clean, GitHub-ready `README.md`** that is tighter, clearer, and better aligned with your actual script (including SLURM, data sources, and full pipeline logic):

---

```markdown
# One-Parent De Novo Mutation Rate Pipeline

Estimate germline **de novo mutation (DNM) rates** using:

- Child **HiFi long reads** (BAM)
- Child **Illumina short reads** (VCF)
- One **parent Illumina VCF**

Designed for scenarios where **only one parent is available**.

---

## Overview

This pipeline performs:

1. Data download (HiFi BAM, Illumina VCF, reference)
2. VCF preprocessing
3. Long-read–guided phasing (Whatshap)
4. Initial DNM detection
5. Long-read validation of candidates
6. Local rephasing around candidate DNMs
7. Callable genome estimation
8. Final mutation rate calculation

---

## Mutation Rate Formula

```

Mutation rate = (# validated DNMs) / (callable genome size)

````

---

## Requirements

### Software

- Python ≥ 3.10
- `bcftools`
- `samtools`
- `whatshap`
- `conda`
- `aws-cli` (for S3 downloads)
- SLURM (recommended)

### Install environment

```bash
conda env create -f environment.yml
conda activate whatshap-env
````

---

## Repository Structure

```
.
├── main.sh                # Main pipeline script
├── environment.yml       # Conda environment
├── data/
│   └── source.txt        # Input data configuration
├── src/                  # Python helper scripts
└── README.md
```

---

## Configuration

Edit:

```
data/source.txt
```

This file defines:

* S3 paths for BAM/VCF
* Reference genome location
* File names

No hardcoded dataset paths inside `main.sh`.

---

## Usage

### Run full pipeline

```bash
bash main.sh \
  --all \
  --prj-dir /path/to/project \
  --sample-child NA12879 \
  --sample-parent NA12878
```

---

### Run specific parts

```bash
bash main.sh \
  --part 1 2 2b \
  --prj-dir /path/to/project \
  --sample-child NA12879 \
  --sample-parent NA12878
```

---

## Pipeline Parts

| Part | Description                     |
| ---- | ------------------------------- |
| 0    | Print configuration             |
| 1    | Data download + preprocessing   |
| 2    | Child phasing (Whatshap)        |
| 2b   | Initial DNM detection           |
| 3    | Local rephasing                 |
| 3b   | Refined DNM detection           |
| 4    | Callable genome + mutation rate |
| 5    | Cleanup                         |

---

## Key Parameters

| Parameter          | Description                 | Default |
| ------------------ | --------------------------- | ------- |
| `--cpus`           | CPUs for phasing            | 8       |
| `--time`           | SLURM walltime              | 6:00:00 |
| `--min-rdepth`     | Min read depth              | 15      |
| `--max-rdepth`     | Max read depth              | 50      |
| `--gt-qual`        | Genotype quality            | 30      |
| `--nv-quantile`    | Variant density threshold   | 0.75    |
| `--mm-diff-min`    | Mismatch threshold          | 0.1     |
| `--window`         | Local rephasing window (bp) | 20000   |
| `--alt-read-count` | Min ALT-supporting reads    | 8       |

---

## Output

All outputs are written to:

```
<PRJ_DIR>/<child>_phasedvcf/
```

### Final results

* `final_dnmc_<child>-from-<parent>.tsv`
* `callable_genome.txt`
* Mutation rate (printed to stdout)

---

## Pipeline Architecture

```
PART 1   → Data download + preprocessing
PART 2   → Long-read phasing
PART 2b  → Initial DNM detection
PART 3   → Local rephasing
PART 3b  → Refinement
PART 4   → Callable genome + mutation rate
PART 5   → Cleanup
```

---

## Notes

* Uses **haplotype block inconsistency** to detect DNMs
* Long reads are used for:

  * Phasing
  * Validation
* Parent is assumed **reference homozygous** at true DNMs
* Numerator and denominator use consistent filtering

---

## Reproducibility

* Modular execution (`--part`)
* Deterministic filtering
* Explicit callable genome estimation
* External data fully configurable via `source.txt`

---



```

