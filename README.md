Perfect — for **public GitHub**, the README should be:

* Clear in 30 seconds
* Installable in 2 minutes
* Runnable in 1 copy-paste block
* Scientifically precise but not verbose
* Architecture-light, execution-heavy

Below is a **clean, public-facing README** tailored for your pipeline.

---

# One-Parent De Novo Mutation Rate Pipeline

Estimate germline **de novo mutation (DNM) rates** using:

* Child **HiFi long reads** (BAM)
* Child **Illumina short reads** (VCF)
* One **parent Illumina short-read** dataset

Designed for scenarios where only **one parent is sequenced**.

---

## Method Overview

1. Phase child short reads using long reads (Whatshap)
2. Merge phased child with unphased parent
3. Detect haplotype mismatches (DNM candidates)
4. Validate candidates with long-read support
5. Rephase locally around candidates
6. Estimate callable genome
7. Compute mutation rate

[
\text{Mutation rate} = \frac{\text{Validated DNMs}}{\text{Callable genome size}}
]

---

## Installation

### 1️⃣ Clone repository

```bash
git clone https://github.com/<your-username>/one-parent-dnm.git
cd one-parent-dnm
```

### 2️⃣ Create environment

```bash
conda env create -f environment.yml
conda activate platinum-dnm
```

### Requirements

* Python ≥ 3.10
* bcftools
* samtools
* whatshap
* SLURM (for cluster execution)
* AWS CLI (if downloading from S3)

---

## Configuration

Edit:

```
data/source.txt
```

This file defines:

* S3 locations
* BAM / VCF filenames
* Reference filename

No hardcoded paths inside `main.sh`.

---

## Running the Pipeline

### Run full workflow (not recommended yet)

```bash
./main.sh \
  --all \
  --prj-dir /path/to/project \
  --sample-child NA12879 \
  --sample-parent NA12878
```

---

### Run specific stages

| Part | Description                     |
| ---- | ------------------------------- |
| 0    | Print configuration             |
| 1    | Download + preprocessing        |
| 2    | Child phasing                   |
| 2b   | Initial DNM detection           |
| 3    | Local rephasing                 |
| 3b   | Refined DNM detection           |
| 4    | Callable genome + mutation rate |
| 5    | Cleanup                         |

Example:

```bash
./main.sh \
  --part 1 2 \
  --prj-dir /path/to/project \
  --sample-child NA12879 \
  --sample-parent NA12878
```

---

## Output

Final results are stored in:

```
PRJ_DIR/<child>_phasedvcf/
```

Key outputs:

* `final_dnmc_<child>-from-<parent>.tsv`
* `callable_genome.txt`
* Mutation rate summary (printed at completion)

---

## Runtime Directory Structure

```
PRJ_DIR/
├── reference/
├── hifi/
└── illumina-dragen/
```

Directories are created automatically if missing.

---

## Scientific Assumptions

* Accurate long-read phasing
* Parent homozygous reference at true DNM sites
* Identical filtering for numerator and denominator
* Phase blocks reflect true haplotypes

---

## Reproducibility

* Modular execution (`--part`)
* Deterministic filtering
* Explicit denominator calculation
* Long-read validation stage
* Dataset configuration isolated in `source.txt`

---

## Citation

If using this pipeline in a publication:

> Nguyen, One-Parent De Novo Mutation Rate Pipeline (GitHub)

---


