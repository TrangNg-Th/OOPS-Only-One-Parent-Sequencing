# OOPS — Only One Parent Sequencing

Estimate germline **de novo mutation (DNM) rates** from a trio where only **one parent** has sequencing data.

---

## The idea in one paragraph

Standard DNM calling needs both parents sequenced. OOPS doesn't. It uses the **child's long reads** to phase the child's genome into two haplotypes — one inherited from each parent — then compares each haplotype, SNP by SNP, against the one sequenced parent. The haplotype inherited from the sequenced parent matches them almost everywhere; the other haplotype (from the absent parent) mismatches a lot, as expected. A site on the *matching* haplotype that still doesn't match the parent is a **de novo mutation candidate** — a variant neither parent transmitted. Long reads from a single parent are enough to phase the child and assign these confidently.

---

## What you need

| Sample | Data |
|---|---|
| Child | Long-read BAM (HiFi **or** ONT) |
| Sequenced parent | Illumina BAM |
| Both (child + sequenced parent) | One **joint** Illumina VCF that already contains both samples' genotypes (`FORMAT/DP`, `GQ`, `GT`, etc.) — Part 1c splits it per-sample. Neither sample needs a raw Illumina BAM beyond what's listed above; genotype-level Illumina evidence comes entirely from this VCF, and the parent's Illumina BAM is used separately for read-level ALT confirmation (Part 2b/6a). |
| Absent parent | — nothing — |
| Reference | The FASTA all BAMs were mapped to |

**Dependencies:** `bcftools`, `samtools`, `tabix`, `whatshap`, `wget`, `pysam`, `pandas`, `numpy`, `scipy`, `matplotlib`, `aws-cli` (S3), SLURM.

```bash
git clone https://github.com/TrangNg-Th/OOPS-Only-One-Parent-Sequencing && cd OOPS-Only-One-Parent-Sequencing

## If you want to build a conda environment from scratch
conda activate build-env
conda-build recipe/ -c conda-forge -c bioconda

```



**Already have an environment named `oops` for something else?** Every phasing/haplotagging step in `main.sh` runs `conda activate "${CONDA_ENV_NAME}"` internally (default `oops`) before calling `whatshap`/`bcftools`/`samtools`. Pass `--conda-env <name>` to point it at whatever you actually named the environment that has these tools installed. For ex: `--conda-env oops`

---

## Pipeline at a glance

```
child long reads + child Illumina VCF
        │  phase child into haplotype blocks (WhatShap)
        ▼
  compare each haplotype to the parent's genotype, per phase block
        │  one haplotype matches the parent, the other doesn't
        ▼
  a single mismatch on the matching haplotype = DNM candidate
        │  confirm ALT allele sits on one haplotype only (long reads)
        ▼
  local rephasing removes phase-switch false positives
        │
        ▼
  callable-genome estimate = the denominator
        │
        ▼
  mutation rate = validated DNMs / callable genome
```

---

## Configuration

All dataset-specific paths live in a **source file** under `data/` — nothing is hardcoded in `main.sh`. Copy the template and edit it:

```bash
cp data/example_readtype-coverage_parent_child.txt \
   data/<datatype>_<coverage>_<child>_<parent>.txt
```

It defines `BAM_CHILD_URL`/`NAME_BAM_CHILD`/`NAME_BAMIDX_CHILD` (child long-read BAM), `BAM_PARENT_URL`/`NAME_BAM_PARENT`/`NAME_BAMIDX_PARENT` (parent Illumina BAM), `VCF_PATH`/`NAME_VCF_FILE` (the one joint Illumina VCF covering both samples), and `REFERENCE_PATH`/`NAME_REFERENCE_FILE`/`NAME_REFERENCE`. Without `--source`, the example file is used.

**Slurm account / partition:** every job `main.sh` submits bills `-A r00379` by default — it's the our HPC allocation but is certainly not for you. Pass `--account <your-account>` (required on most clusters) and `--partition <your-partition>` (this is optional; omit to let Slurm pick the cluster's default) on every invocation.

---

## Quick start

- Run `--part 0` first to check yoour project's path, and all the parameters involved. 

- Part 1 downloads the child long-read BAM + index, the parent's Illumina BAM + index, the one joint Illumina VCF (genotypes for *both* samples), and the reference (then `faidx`).

- There are two ways to run the pipeline:
      - If you have **not phased** the child using long reads, you can run part 2 to phase the child via Whatshap.
      - If you have **already phased** the child, you can skip part 2 and run part 2b to detect DNM candidates. Make sure that the phased VCF has `FORMAT/PS` and phased genotypes (`|` instead of `/`) for the child sample.

            - Once you have the **phased VCF**, you can run part 2b to detect DNM candidates, then part 3 to locally rephase around candidates, then part 3b to refine DNM detection, and finally part 4 to estimate callable genome and mutation rate.
            
            - Part 2b ->  3-> 3b-> 4 can be run **chained** (recommended) or **one part at a time**. Both ways produce the same final outputs.

**Chained (recommended):** one command submits all four as a Slurm dependency chain — each step waits for the previous to finish.

```bash
bash main.sh --part chain \
  --prj-dir /path/OOPS_hifi_5x_NA12879_NA12878 \
  --sample-child NA12879 --sample-parent NA12878 \
  --source hifi_5x_NA12879_NA12878.txt \
  --nv-quantile 0.5 --alt-read-count 2 --window 20000 --total-rd-ct-min 10 \
  > logs/hifi_5x_NA12879_NA12878.log 2>&1
```

Chain artifacts and per-step logs land in `<prj-dir>/<child>_phasedvcf/chain_2b_4/logs/`. Monitor with `squeue -u $USER`. (If a step fails, `afterok` cancels the downstream steps — that's expected.)


**One part at a time:** same flags, just change `--part` (`0 → 1a → 1b → 1c → 2 → 2b → 3 → 3b → 4`). Every part still runs standalone.

---

## Parts reference

Run order: `0 → 1a → 1b → 1c → 2 → 2b → 3 → 3b → 4`. Parts 5 and 6a are optional refinement.

### 0 — Configuration check
Prints the resolved parameters and dependencies; submits nothing. Verify the child BAM, VCF, and reference paths before continuing.

### 1a / 1b / 1c — Data prep - optional
- **1a** downloads the child long-read BAM + index, the parent's Illumina BAM + index, the one joint Illumina VCF (genotypes for *both* samples), and the reference (then `faidx`).
- **1b** strips any existing `HP` tags from the child's long-read BAM, sorts, indexes → `<child>_clean.bam`. (The parent's Illumina BAM is left as downloaded — it's only read directly, never modified.)
- **1c** splits the joint VCF into the two per-sample VCFs it already contains, removes `PS` tags, indexes → per-sample `*.unphased.noPS.vcf.gz`. No separate child Illumina BAM is downloaded or needed — child genotype/depth/quality data comes entirely from this VCF.

Each submits a short Slurm job.

### 2 — Initial phasing
Submits `build_hapl_<child>`: filters the child VCF to diploid sites, runs `whatshap phase` (producing PS-tagged blocks), then `whatshap haplotag` to stamp every read with `HP1`/`HP2` → `<child>_HP.bam`.
*Check `*.stats.tsv`: at 50× HiFi expect block N50 >100 kb and >85% of hets phased.*

**If you're using a different phasing program:** pass `--external-phased-vcf <FILE>` to skip `whatshap phase` and use a VCF already phased by another tool instead (**WARNING:** your vcf file must carry `FORMAT/PS` and phased `|` genotypes). However, `whatshap` itself is not optional — only the initial phase call is swappable. Haplotagging and phase-block stats always run via `whatshap` afterward regardless of this flag, since `whatshap haplotag`/`stats` accept any correctly phased VCF no matter which tool produced it. 

Either way, Part 2 — and again Part 2b, independently — checks the phased VCF for a non-missing `FORMAT/PS` and fails fast with a clear error if phasing didn't actually produce phase sets, instead of silently yielding zero candidates deep in Part 2b.

### 2b — First DNM detection (if you want to chain Parts 2b → 3 → 3b → 4, skip this and run `--part chain` instead)
Merges phased child + unphased parent VCFs, keeps PS-tagged child SNPs, counts H0/H1 mismatches per block, flags blocks where exactly one haplotype has a single mismatch, and validates each candidate in the haplotagged long reads (ALT must sit on one haplotype, ≥ `--alt-read-count`, ≥ `--total-rd-ct-min` total).
Outputs in `<child>_phasedvcf/mismatch_analysis/`: `*_mismatch.tsv`, `*_dnmc.bed/.tsv`, `<child>_LR_validated_dnmc.bed`, `<child>_dnmc_plusminus<W>kb.bed`.

**This part is included in the `--part chain` flow**; you only need to run it standalone if you want to inspect the initial candidates before local rephasing.

### 3 — Local rephasing around candidates
Submits a job that re-runs WhatShap in a ±`--window` bp region (default 20kb) around each candidate. Global phasing accumulates switch errors over distance; local rephasing gives a cleaner haplotype call near each site. A candidate that vanishes here was a phase-switch artifact.

This step always rephases with `whatshap phase` — it's not swappable via a flag. If you'd rather phase these windows with your own program, let Part 3 run once to lay out its inputs, then before running Part 3b (or `--part chain`), replace its output file with your own phased VCF for the same regions (must carry `FORMAT/PS` and phased genotypes):

- Region BED Part 3 outputs: `<prj-dir>/<child>_phasedvcf/mismatch_analysis/<child>_dnmc_plusminus<window-kb>kb.bed`

**This part is included in the `--part chain` flow**; you only need to run it standalone if you want to inspect the initial candidates before local rephasing.

### 3b — Refined DNM detection
Re-merges the locally rephased child VCF with the parent, re-extracts phased SNPs in the candidate windows, and re-counts with the same asymmetry logic — a confirmation pass. Survivors are promoted.
Output: `<child>_phasedvcf/final_dnmc_<child>-from-<parent>.tsv`.

**This part is included in the `--part chain` flow**; you only need to run it standalone if you want to inspect the initial candidates before local rephasing.

### 4 — Callable genome + mutation rate
Estimates the denominator instead of using raw genome size: samples het SNPs from the same blocks, applies the **same** depth/GQ/long-read filters used for DNMs, and extrapolates.

**This part is included in the `--part chain` flow**; you only need to run it standalone if you want to inspect the initial candidates before local rephasing.

```
callable_genome = (SNPs_qualified / SNPs_sampled) × accessible_bases
mutation_rate   = N_validated_DNMs / callable_genome
```

Prints the summary and writes `mismatch_analysis/denum_calcul/callable_genome.txt`.

### 5 — High-mismatch block rephasing (optional)
Use if the rate looks inflated or candidates cluster in one block. Finds blocks where **both** haplotypes mismatch heavily (symmetric = switch error, not signal), splits them into `--window_rephase` sub-windows, and rephases each via a Slurm array (≤200 tasks) + a summary job.

### 6a — Re-evaluate DNMs in rephased blocks (optional)
Re-runs the full detection logic on Part 5's rephased blocks and computes an independent callable genome for them, then reports a candidate-weighted rate combining the original and rephased sets. Run after Part 5 finishes; `--part 6a` auto-submits its internal merge + dependent steps.

---

## Key parameters

| Parameter | Controls | Default | Change when |
|---|---|---|---|
| `--min-rdepth` / `--max-rdepth` | Depth window per site | 15 / 50 | Shallow data → lower min; high cov → raise max |
| `--gt-qual` | Min genotype quality | 30 | Many sites dropped unexpectedly |
| `--nv-quantile` | Min block size (SNP-count quantile) | 0.5 | Raise to 0.75 to cut short-block false positives |
| `--mm-diff-min` | Min H0/H1 mismatch asymmetry | 0.1 | Raise in noisy data |
| `--alt-read-count` | Min ALT-supporting long reads | 8 | Lower to 1–3 for ≤10× coverage |
| `--total-rd-ct-min` | Min total (REF+ALT) reads | 5 | Raise for higher coverage |
| `--window` | Region for local rephasing (bp) | 20000 | Raise for large phase blocks |
| `--window_rephase` | Sub-window in Part 5 (bp) | 100000 | Lower for very long blocks |
| `--threshold_rephase` | Mismatch count to flag a block | 100 | Lower to rephase more aggressively |
| `--cpus` / `--time` | Phasing job resources | 8 / 10:00:00 | Large BAMs / many candidates |
| `--account` | Slurm account to bill (`-A`) | `r00379` (the author's — change this) | Always, unless you happen to share the author's allocation |
| `--partition` | Slurm partition/queue | unset (cluster default) | Your cluster requires an explicit partition |
| `--conda-env` | Env name `conda activate`d before whatshap/bcftools/samtools calls | `oops` | You already have an env named `oops` for something else |

---

## Main outputs

| File | Description |
|---|---|
| `final_dnmc_<child>-from-<parent>.tsv` | Final validated DNMs |
| `mismatch_analysis/denum_calcul/callable_genome.txt` | Callable genome + denominator breakdown |
| `mismatch_analysis/<child>_LR_validated_dnmc.bed` | Long-read-validated candidates (initial pass) |
| `mismatch_analysis/<parent>_<child>_mismatch.tsv` | Per-block H0/H1 mismatch counts |
| `<child>.illumVCF_LRbam.blocks.tsv` | Phase-block coordinates and SNP counts |
| `rephased_blocks/rephase_summary.txt` | Part 5 array success/failure summary |

---

## Repository layout

```
main.sh                  # pipeline entry point
environment.yml          # conda env
data/                    # source config templates
src/
  count_mismatches.py            # per-block mismatch counting + candidate ID
  count_rephase_mismatches.py    # same logic on rephased blocks (6a)
  dnmc_readcheck.py              # long-read validation (ALT on one haplotype)
  parent_readcheck.py            # validates candidates against parent Illumina BAM
  callable_genome.py             # SNP sampling + callable-genome estimate (4)
  rephase_blocks.py              # flags high-mismatch blocks for Part 5
  fix_PhaseSet.py                # reformats WhatShap PS output to TSV
  dnm_vs_depth.py                # standalone plot: DNM calls/rate vs read depth
  slurm_scripts_helper/          # optional one-off SLURM utilities (add read groups, check/downsample BAM coverage)
```

---

## Limitations

- **Germline only** — somatic mutations are out of scope.
- The sequenced parent should be **homozygous reference** at true DNM sites; heterozygous parent sites are filtered but reduce sensitivity.
- **Phase-switch errors** are the main false-positive source; Parts 5/6a mitigate but don't eliminate them.
- At low long-read coverage, some candidates are genuinely unresolvable per site 
- Validated on HiFi, ONT from one trio (NA12879/NA12878/NA12877).Multi-pedigree validation are ongoing.


---
### TODOs
- Update pipeline to reproduce work for PacBio Hifi on trio NA12878/NA12879/NA12877
- 