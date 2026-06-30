#!/bin/bash

set -eo pipefail
export PS1=${PS1:-}  


################################################################################
# Project: OOPS - Only One Parent Sequencing Mutation Call Pipeline
################################################################################


usage() {
cat << 'EOF'

Only One Parent Sequencing (OOPS) — Mutation Call Pipeline
==========================================================

Estimate germline de novo mutation (DNM) rates from a trio in which only one
parent has long-read data.

------------------------------------------------------------------
USAGE
------------------------------------------------------------------
  ./main.sh --part <PART> --prj-dir <DIR> --sample-child <ID> \\
            --sample-parent <ID> --source <FILE> [options]

------------------------------------------------------------------
REQUIRED
------------------------------------------------------------------
  --prj-dir <DIR>         Project root directory
  --sample-child <ID>     Child sample ID  (e.g. NA12879)
  --sample-parent <ID>    Parent sample ID (e.g. NA12878)

------------------------------------------------------------------
WHAT TO RUN  (--part)
------------------------------------------------------------------
  Normal order:  0 -> 1a -> 1b -> 1c -> 2 -> 2b -> 3 -> 3b -> 4
  Optional:      5 -> 6a   (phase-switch refinement)

  0      Print resolved configuration and exit (run this first)
  1a     Download data (child long-read BAM, parent Illumina BAM,
         joint Illumina VCF, reference)
  1b     BAM preprocessing (strip HP tags, sort, index)
  1c     VCF preprocessing (split joint VCF by sample, drop PS, index)
  2      Initial phasing (WhatShap) + haplotag the child BAM
  2b     First DNM detection
  3      Local rephasing around candidate DNMs (submits a job)
  3b     Refined DNM detection
  4      Callable-genome estimation + final mutation rate
  5      Flag high-mismatch blocks and locally rephase them (optional)
  6a     Re-evaluate DNMs in rephased blocks (optional)

  chain  Submit 2b -> 3 -> 3b -> 4 as ONE dependency chain of Slurm jobs
         (each step waits for the previous to finish). Individual parts can
         still be run by hand exactly as above.

  --all  Run: 0 1a 1b 1c 2 2b 3 3b 4 5 6a

  --part accepts several tokens:  --part 1a 1b 1c 2

------------------------------------------------------------------
OPTIONS  (defaults in brackets)
------------------------------------------------------------------
  Resources
    --cpus <N>              CPUs for phasing jobs           [8]
    --time <HH:MM:SS>       Slurm wall time                 [10:00:00]
    --reference <FILE>      Reference FASTA                 [<prj-dir>/reference/chm13v2.0_maskedY_rCRS.fa]
    --source <FILE>         Source config under data/       [example_readtype-coverage_parent_child.txt]
    --account <NAME>        Slurm account to bill (-A)      [r00379]
                             Default is the original author's HPC allocation
                             -- almost certainly wrong for you. Set this to
                             your own Slurm account before running anything.
    --partition <NAME>      Slurm partition/queue           [unset]
                             Unset lets Slurm pick the cluster's default
                             partition. Set this if your cluster requires an
                             explicit partition/queue name.
    --conda-env <NAME>      Conda env with whatshap/bcftools/samtools  [oops]
                             Every phasing/haplotagging/stats step runs
                             `conda activate <NAME>` before calling whatshap.
                             Set this to whatever you named the environment
                             when you ran `conda env create -f environment.yml`
                             (or created when installing the conda-build
                             package) -- e.g. pass `--conda-env whatshap-env`
                             if you already have an unrelated environment
                             named "oops" and put this pipeline's tools under
                             a different name instead.

  Phasing (Part 2)
    --external-phased-vcf <FILE>  Skip `whatshap phase`; use a VCF already
                                   phased by another tool instead. The file
                                   MUST carry FORMAT/PS and phased (|)
                                   genotypes -- Part 2 checks for this and
                                   fails fast if it doesn't. Haplotagging and
                                   phase-block stats ALWAYS run via whatshap
                                   afterward, regardless of this flag --
                                   whatshap is not an optional dependency,
                                   only the initial phase call is swappable.
                                   [unset -> run whatshap phase as before]

  Local rephasing (Part 3)
    Part 3 always rephases candidate-window regions with `whatshap phase`
    (this refinement step is not swappable via a flag). If you'd rather
    phase these windows with your own program, run Part 3 to let it lay out
    the inputs, then before running Part 3b/chain replace its output file:
      <prj-dir>/<child>_phasedvcf/mismatch_analysis/<child>.illumVCF_LRbam.phased.<window-kb>kb.vcf.gz
    with your own phased VCF for that same region set (must carry FORMAT/PS
    and phased genotypes for <child>). The region BED Part 3 phases is:
      <prj-dir>/<child>_phasedvcf/mismatch_analysis/<child>_dnmc_plusminus<window-kb>kb.bed

  Site / genotype filters
    --min-rdepth <N>        Min read depth                  [15]
    --max-rdepth <N>        Max read depth                  [50]
    --gt-qual <N>           Min genotype quality (GQ)       [30]

  Block selection
    --nv-quantile <F>       Min block size by SNP-count quantile  [0.5]
    --mm-diff-min <F>       Min H0/H1 mismatch asymmetry          [0.1]

  Long-read validation
    --min-base-qual <N>     Min base quality                [20]
    --min-map-qual <N>      Min mapping quality             [20]
    --window <BP>           Window around each candidate    [20000]
    --alt-read-count <N>    Min ALT-supporting reads        [8]  (use 1-3 for <=10x)
    --total-rd-ct-min <N>   Min total (REF+ALT) reads       [5]
    --verbose <T|F>         Verbose LR validation           [T]

  Phase-switch refinement (Parts 5/6a)
    --window_rephase <BP>     Sub-window when splitting blocks  [100000]
    --threshold_rephase <N>   Min mismatch count to flag block  [100]

  Misc
    --exclude-chroms <LIST> Comma-separated chroms to skip  [chrX,chrY,chrM]

------------------------------------------------------------------
EXAMPLES
------------------------------------------------------------------
  # 0. Check configuration
  ./main.sh --part 0 --prj-dir /path/OOPS_hifi_5x_NA12879_NA12878 \\
    --sample-child NA12879 --sample-parent NA12878 \\
    --source hifi_5x_NA12879_NA12878.txt

  # Run detection-to-rate as one chained submission
  ./main.sh --part chain --prj-dir /path/OOPS_hifi_5x_NA12879_NA12878 \\
    --sample-child NA12879 --sample-parent NA12878 \\
    --source hifi_5x_NA12879_NA12878.txt \\
    --nv-quantile 0.5 --alt-read-count 2 --window 20000 --total-rd-ct-min 10

  # Or one part at a time
  ./main.sh --part 2b --prj-dir /path --sample-child NA12879 \\
    --sample-parent NA12878 --source hifi_5x_NA12879_NA12878.txt

EOF
exit 1
}

##############################################
# Default values
##############################################

CPUS=8
TIME="10:00:00"
MIN_RDEPTH=15
MAX_RDEPTH=50
GT_QUAL=30
NV_QUANTILE=0.5
MM_DIFF_MIN=0.1
MIN_BASE_QUAL=20
MIN_MAP_QUAL=20
WINDW=20000
ALT_READ_COUNT=8
VERBOSE=T
RECOUNT=T
NOTRECOUNT=F
THRESHOLD_REPHASE=100
EXCLUDE_CHROMS="chrX,chrY,chrM"
WINDOW_REPHASE=100000
REPHASE_SUFFIX="_rephase"
TOTAL_READ_COUNT_MIN=5
EXTERNAL_PHASED_VCF=""
ACCOUNT="r00379"
PARTITION=""
CONDA_ENV_NAME="oops"


##############################################
# Parse arguments
##############################################

if [[ $# -eq 0 ]]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        --part)
            shift
            while [[ $# -gt 0 && $1 != --* ]]; do
                PARTS+=("$1")
                shift
            done
            ;;
        --all)
            PARTS=("0" "1a" "1b" "1c" "2" "2b" "3" "3b" "4" "5" "6a" "chain")
            shift
            ;;

        --prj-dir) PRJ_DIR="$2"; shift 2 ;;
        --reference) REF="$2"; shift 2 ;;
        --sample-child) SAMPLE_CHILD="$2"; shift 2 ;;
        --sample-parent) SAMPLE_PARENT="$2"; shift 2 ;;
        --cpus) CPUS="$2"; shift 2 ;;
        --time) TIME="$2"; shift 2 ;;
        --min-rdepth) MIN_RDEPTH="$2"; shift 2 ;;
        --max-rdepth) MAX_RDEPTH="$2"; shift 2 ;;
        --gt-qual) GT_QUAL="$2"; shift 2 ;;
        --nv-quantile) NV_QUANTILE="$2"; shift 2 ;;
        --mm-diff-min) MM_DIFF_MIN="$2"; shift 2 ;;
        --min-base-qual) MIN_BASE_QUAL="$2"; shift 2 ;;
        --min-map-qual) MIN_MAP_QUAL="$2"; shift 2 ;;
        --window) WINDW="$2"; shift 2 ;;
        --alt-read-count) ALT_READ_COUNT="$2"; shift 2 ;;
        --verbose) VERBOSE="$2"; shift 2 ;;
        --source) SOURCE="$2"; shift 2 ;;
        --window_rephase) WINDOW_REPHASE="$2"; shift 2 ;;
        --threshold_rephase) THRESHOLD_REPHASE="$2"; shift 2 ;;
        --exclude-chroms) EXCLUDE_CHROMS="$2"; shift 2 ;;
        --total-rd-ct-min) TOTAL_READ_COUNT_MIN="$2"; shift 2 ;;
        --external-phased-vcf) EXTERNAL_PHASED_VCF="$2"; shift 2 ;;
        --account) ACCOUNT="$2"; shift 2 ;;
        --partition) PARTITION="$2"; shift 2 ;;
        --conda-env) CONDA_ENV_NAME="$2"; shift 2 ;;
        --help) usage ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
done

##############################################
# Required parameter check
##############################################

if [[ -z "${PRJ_DIR:-}" || -z "${SAMPLE_CHILD:-}" || -z "${SAMPLE_PARENT:-}" ]]; then
    echo "ERROR: Missing required arguments."
    usage
fi



##############################################
# Derived variables
##############################################

ILLUM_DIR=${PRJ_DIR}/illumina-dragen
BAM_DIR=${PRJ_DIR}/bam
REF_DIR=${PRJ_DIR}/reference
v=$((WINDW / 1000))



# Make directory
mkdir -p ${PRJ_DIR}
mkdir -p ${ILLUM_DIR}
mkdir -p ${BAM_DIR}
mkdir -p ${REF_DIR}




##############################################
# Load data source configuration
##############################################
# source.txt lives in data/ same level as this script main.sh
WORKING_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEPENDENCIES_FILE="${WORKING_DIR}/environment.yml"


##############################################
# Check if the argument --source is given, if not, use default source file
##############################################
if [[ -n "${SOURCE:-}" ]]; then
    SOURCE_FILE="${WORKING_DIR}/data/${SOURCE}"
else
    SOURCE_FILE="${WORKING_DIR}/data/example_readtype-coverage_parent_child.txt"
fi



## Read the source file
if [[ ! -f "${SOURCE_FILE}" ]]; then
    echo "ERROR: ${SOURCE} not found at ${SOURCE_FILE}"
    exit 1
fi

set -a
source <(grep -v '^##' "${SOURCE_FILE}")
set +a

####################################################
## Final BAM / VCF paths
####################################################
BAM_CHILD="${BAM_CHILD_URL}/${NAME_BAM_CHILD}"
BAMIDX_CHILD="${BAM_CHILD_URL}/${NAME_BAMIDX_CHILD}"

# VCF
VCF_FULL_PATH="${VCF_PATH}/${NAME_VCF_FILE}"
NAME_VCF="${ILLUM_DIR}/${NAME_VCF_FILE}"

# Reference
REF="${REF:-${REF_DIR}/${NAME_REFERENCE_FILE}}"

# Long-read BAM naming
ORIG_BAM_PATH="${BAM_DIR}/${NAME_BAM_CHILD}"
ORIG_BAM_BASENAME="$(basename "${NAME_BAM_CHILD}" .bam)"
CLEAN_BAM_PATH="${BAM_DIR}/${ORIG_BAM_BASENAME}_clean.bam"
CLEAN_BAM_INDEX="${CLEAN_BAM_PATH}.bai"
HP_BAM_PATH="${BAM_DIR}/${ORIG_BAM_BASENAME}_HP.bam"
HP_BAM_INDEX="${HP_BAM_PATH}.bai"

# Parent Illumina BAM (for parental allele evidence at DNM candidates)
ORIG_BAM_PARENT_PATH="${BAM_DIR}/${NAME_BAM_PARENT}"

# If the user wants to set partition, if not, blank 
SBATCH_PARTITION_LINE=""
if [[ -n "${PARTITION}" ]]; then
    SBATCH_PARTITION_LINE="#SBATCH --partition=${PARTITION}"
fi

##############################################
# Print configs
##############################################

echo "================ Pipeline Configuration ================"
echo "Project directory                     : ${PRJ_DIR}"
echo "Child sample                          : ${SAMPLE_CHILD}"
echo "Parent sample                         : ${SAMPLE_PARENT}"
echo "CPUs                                  : ${CPUS}"
echo "Walltime                              : ${TIME}"
echo "Min read depth (SR)                   : ${MIN_RDEPTH}"
echo "Max read depth (SR)                   : ${MAX_RDEPTH}"
echo "Genotype quality cutoff(SR)           : ${GT_QUAL}"
echo "Quantile of nbr variants per blocks   : ${NV_QUANTILE}"
echo "Mismatch diff threshold               : ${MM_DIFF_MIN}"
echo "Min base quality (LR)                 : ${MIN_BASE_QUAL}"
echo "Min mapping quality (LR)              : ${MIN_MAP_QUAL}"
echo "Window size (bp)                      : ${WINDW}"
echo "Window size for rephasing (bp)        : ${WINDOW_REPHASE}"
echo "Mismatch count threshold for rephase  : ${THRESHOLD_REPHASE}"
echo "Alt read count (LR)                   : ${ALT_READ_COUNT}"
echo "Verbose LR validation                 : ${VERBOSE}"
echo "Data source file name                 : ${SOURCE_FILE}"
echo "Slurm account (-A)                    : ${ACCOUNT}"
echo "Slurm partition                       : ${PARTITION:-<cluster default>}"
echo "Conda environment (phasing/whatshap)  : ${CONDA_ENV_NAME}"
echo "Bam file of the parent                : ${BAM_PARENT_URL}/${NAME_BAM_PARENT}"
echo "Bam file of the child                 : ${BAM_DIR}/${NAME_BAM_CHILD}"
echo "Vcf file                              : ${NAME_VCF}"
echo "Reference file                        : ${REF}"
echo "Excluded chroms                       : ${EXCLUDE_CHROMS}"
echo "Total read count minimum for dnmc     : ${TOTAL_READ_COUNT_MIN}"
echo "Parts to run                          : ${PARTS[*]}"
echo "END SUMMARY"
echo "========================================================="
echo ""


##############################################
# Load modules
##############################################
# If on HPS, load these modules, if not, skip

# ===========================================================================
## STEP 0 - Load modules and print configuration summary
## ==========================================================================
load_modules() {

    if command -v module &> /dev/null; then
        module load bcftools || true
        module load samtools || true
        module load aws-cli || true
        module load conda || true
    fi

}

# ===========================================================================
## STEP 1A - Download data (BAM, VCF, reference)
## ==========================================================================

download_data_job() {

    local OUT=./download_data_${SAMPLE_CHILD}.slurm
    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J download_data_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH -A ${ACCOUNT}
${SBATCH_PARTITION_LINE}
#SBATCH -o slurm_output/download_.%j.txt
#SBATCH -e slurm_output/download_.%j_.err

set -euo pipefail

if command -v module &> /dev/null; then
    module load aws-cli/2.25.5 || true
    module load samtools || true
fi

SAMPLE=${SAMPLE_CHILD}
PRJ_DIR=${PRJ_DIR}
BAM_DIR=${BAM_DIR}
ILLUM_DIR=${ILLUM_DIR}


mkdir -p \${BAM_DIR}
mkdir -p \${ILLUM_DIR}


echo "[download] Downloading long-read BAM for \${SAMPLE}"

aws s3 cp --no-sign-request ${BAM_CHILD} \${BAM_DIR}/${NAME_BAM_CHILD}

aws s3 cp --no-sign-request ${BAMIDX_CHILD} \${BAM_DIR}/${NAME_BAMIDX_CHILD}

echo "[download] Downloading parent Illumina BAM for ${SAMPLE_PARENT}"

aws s3 cp --no-sign-request ${BAM_PARENT_URL}/${NAME_BAM_PARENT} \${BAM_DIR}/${NAME_BAM_PARENT}
aws s3 cp --no-sign-request ${BAM_PARENT_URL}/${NAME_BAMIDX_PARENT} \${BAM_DIR}/${NAME_BAMIDX_PARENT}

if [[ ! -f "\${BAM_DIR}/${NAME_BAM_PARENT}.bai" \
   && ! -f "\${BAM_DIR}/$(basename ${NAME_BAM_PARENT} .bam).bai" ]]; then
    echo "[download] Parent BAM index not found; building one"
    samtools index \${BAM_DIR}/${NAME_BAM_PARENT}
fi

mkdir -p "${REF_DIR}"

echo "[download] Downloading Reference file"

wget -P "${REF_DIR}" "${REFERENCE_PATH}"

gunzip "${REF_DIR}"/"${NAME_REFERENCE_FILE}"

samtools faidx "${REF_DIR}"/"${NAME_REFERENCE_FILE}"

echo "[download] Downloading the Illumina VCF"

aws s3 cp --no-sign-request "${VCF_FULL_PATH}" "${NAME_VCF}"

echo "[download] Done"
EOF

    chmod +x "${OUT}"

    ## Extract the id of the job

    JOBID=$(sbatch --parsable "${OUT}")
    echo "${JOBID}"
    rm "${OUT}"
}


# ===========================================================================
## STEP 1B - BAM preprocessing (remove HP tags if present in bam, sort, index)
## ==========================================================================

generate_bam_preprocessing_job() {

    local OUT=./preprocess_bam_${SAMPLE_CHILD}.slurm

    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J preprocess_bam_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH -A ${ACCOUNT}
${SBATCH_PARTITION_LINE}
#SBATCH -o slurm_output/prepr_bam_.%j.txt
#SBATCH -e slurm_output/prepr_bam_.%j_.err

set -euo pipefail
if command -v module &> /dev/null; then
    module load samtools || true
fi

INPUT_BAM="${ORIG_BAM_PATH}"
OUTPUT_BAM="${CLEAN_BAM_PATH}"

echo "[bam-preprocess] Input BAM: \${INPUT_BAM}"
echo "[bam-preprocess] Output clean BAM: \${OUTPUT_BAM}"

if [[ ! -f "\${INPUT_BAM}" ]]; then
    echo "ERROR: input BAM not found: \${INPUT_BAM}"
    exit 1
fi

echo "[bam-preprocess] Checking whether BAM contains HP tags"

HAS_HP=0
if samtools view "\${INPUT_BAM}" \
    | awk 'index(\$0, "\tHP:i:") {found=1} END{exit !found}'; then
    HAS_HP=1
fi

if [[ "\${HAS_HP}" -eq 1 ]]; then
    echo "[bam-preprocess] HP tag detected. Removing HP tags from BAM."

    samtools view -h "\${INPUT_BAM}" \
      | awk 'BEGIN{OFS="\t"}
             /^@/ {print; next}
             {
               out = \$1
               for (i = 2; i <= NF; i++) {
                 if (\$i !~ /^HP:i:/) out = out OFS \$i
               }
               print out
             }' \
      | samtools view -b -o "\${OUTPUT_BAM}" -

else
    echo "[bam-preprocess] No HP tag detected. Renaming BAM to clean BAM."
    mv "\${INPUT_BAM}" "\${OUTPUT_BAM}"

    if [[ -f "\${INPUT_BAM}.bai" ]]; then
        mv "\${INPUT_BAM}.bai" "\${OUTPUT_BAM}.bai" || true
    fi
fi

echo "[bam-preprocess] Sorting clean BAM"
samtools sort -o "\${OUTPUT_BAM}.sorted" "\${OUTPUT_BAM}"
mv "\${OUTPUT_BAM}.sorted" "\${OUTPUT_BAM}"

echo "[bam-preprocess] Indexing clean BAM"
samtools index "\${OUTPUT_BAM}"

echo "[bam-preprocess] Done"
EOF

    chmod +x "${OUT}"
    JOBID=$(sbatch --parsable "${OUT}")
    echo "${JOBID}"
    rm "${OUT}"
}


# ===========================================================================
## STEP 1C - VCF preprocessing (split joint VCF by sample, remove PS tags, index)
## ==========================================================================

generate_vcf_preprocessing_job() {

    local OUT=./preprocess_vcf_${SAMPLE_CHILD}.slurm
    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J preprocess_vcf_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH -A ${ACCOUNT}
${SBATCH_PARTITION_LINE}
#SBATCH -o slurm_output/preprocess_vcf_.%j.txt
#SBATCH -e slurm_output/preprocess_vcf_.%j_.err


set -euo pipefail
export PS1=\${PS1:-}

if command -v module &> /dev/null; then
    module load bcftools || true
fi

# Split VCFs
bcftools view -s ${SAMPLE_PARENT} -Oz -o ${ILLUM_DIR}/${SAMPLE_PARENT}.${NAME_REFERENCE}.illumina.vcf.gz ${ILLUM_DIR}/${NAME_VCF_FILE}
bcftools view -s ${SAMPLE_CHILD} -Oz -o ${ILLUM_DIR}/${SAMPLE_CHILD}.${NAME_REFERENCE}.illumina.vcf.gz ${ILLUM_DIR}/${NAME_VCF_FILE}

# Index
tabix -p vcf ${ILLUM_DIR}/${SAMPLE_PARENT}.${NAME_REFERENCE}.illumina.vcf.gz
tabix -p vcf ${ILLUM_DIR}/${SAMPLE_CHILD}.${NAME_REFERENCE}.illumina.vcf.gz

# Clean PS tags (CHILD)
bcftools +setGT ${ILLUM_DIR}/${SAMPLE_CHILD}.${NAME_REFERENCE}.illumina.vcf.gz -Ou -- -t a -n u \
| bcftools annotate -x FORMAT/PS \
| bcftools view -Oz -o ${ILLUM_DIR}/${SAMPLE_CHILD}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz

bcftools index ${ILLUM_DIR}/${SAMPLE_CHILD}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz


# Clean PS tags (PARENT)
bcftools +setGT ${ILLUM_DIR}/${SAMPLE_PARENT}.${NAME_REFERENCE}.illumina.vcf.gz -Ou -- -t a -n u \
| bcftools annotate -x FORMAT/PS \
| bcftools view -Oz -o ${ILLUM_DIR}/${SAMPLE_PARENT}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz

bcftools index ${ILLUM_DIR}/${SAMPLE_PARENT}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz

echo "[vcf-preprocess] Done"

EOF

    chmod +x "${OUT}"
    JOBID=$(sbatch --parsable "${OUT}")
    echo "${JOBID}"
    rm "${OUT}"
}




# ===========================================================================
## STEP 2 - Initial phasing with Whatshap + haplotagging (take about 1 hour per sample)
## ==========================================================================

generate_phasing_job() {
    local OUT="${PRJ_DIR}/build_hapl_${SAMPLE_CHILD}.slurm"

    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J build_hapl_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME}
#SBATCH --mem=12G
#SBATCH -A ${ACCOUNT}
${SBATCH_PARTITION_LINE}
#SBATCH -o slurm_output/build_hapl_%j.txt
#SBATCH -e slurm_output/build_hapl_%j.err

set -euo pipefail
export PS1=\${PS1:-}

if command -v module &> /dev/null; then
    module load bcftools || true
    module load samtools || true
    module load conda || true
fi
conda activate "${CONDA_ENV_NAME}"

REF="${REF}"
SAMPLE="${SAMPLE_CHILD}"

ILLUM_VCF="${ILLUM_DIR}/\${SAMPLE}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz"
CLEAN_BAM="${CLEAN_BAM_PATH}"
HP_BAM="${HP_BAM_PATH}"
OUT_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"

DIPLOID_VCF="\${OUT_DIR}/\${SAMPLE}.illumVCF_LRbam.diploid.vcf.gz"
PHASED_VCF="\${OUT_DIR}/\${SAMPLE}.illumVCF_LRbam.phased.vcf.gz"
PHASE_STATS_TSV="\${OUT_DIR}/\${SAMPLE}.illumVCF_LRbam.stats.tsv"
PHASE_BLOCKS_TSV="\${OUT_DIR}/\${SAMPLE}.illumVCF_LRbam.blocks.tsv"

mkdir -p "\${OUT_DIR}"
mkdir -p slurm_output

echo "[phase] Input VCF         : \${ILLUM_VCF}"
echo "[phase] Clean BAM         : \${CLEAN_BAM}"
echo "[phase] Diploid-only VCF  : \${DIPLOID_VCF}"
echo "[phase] Phased VCF        : \${PHASED_VCF}"
echo "[phase] HP-tagged BAM     : \${HP_BAM}"

if [[ ! -f "\${ILLUM_VCF}" ]]; then
    echo "ERROR: Illumina VCF not found: \${ILLUM_VCF}"
    exit 1
fi

if [[ ! -f "\${CLEAN_BAM}" ]]; then
    echo "ERROR: Clean BAM not found: \${CLEAN_BAM}"
    exit 1
fi

if [[ ! -f "\${CLEAN_BAM}.bai" ]]; then
    echo "[phase] BAM index not found. Creating index for clean BAM."
    samtools index "\${CLEAN_BAM}"
fi

echo "[phase] Removing stale outputs"
rm -f "\${DIPLOID_VCF}" "\${DIPLOID_VCF}.csi" "\${DIPLOID_VCF}.tbi"
rm -f "\${PHASED_VCF}" "\${PHASED_VCF}.csi" "\${PHASED_VCF}.tbi"
rm -f "\${PHASE_STATS_TSV}" "\${PHASE_BLOCKS_TSV}"
rm -f "\${HP_BAM}" "\${HP_BAM}.bai"

EXTERNAL_PHASED_VCF="${EXTERNAL_PHASED_VCF}"

if [[ -z "\${EXTERNAL_PHASED_VCF}" ]]; then

    echo "[phase] Building diploid-only VCF to avoid mixed-ploidy crashes"
    bcftools view -s "\${SAMPLE}" -Ou "\${ILLUM_VCF}" \
      | bcftools view -i 'GT!="hap" && GT!="mis"' \
      -Oz -o "\${DIPLOID_VCF}"

    echo "[phase] Indexing diploid-only VCF"
    bcftools index -f "\${DIPLOID_VCF}"

    echo "[phase] Checking whether any records remain after diploid filtering"
    N_DIPLOID=\$(bcftools index -n "\${DIPLOID_VCF}")
    echo "[phase] Number of records retained: \${N_DIPLOID}"

    if [[ "\${N_DIPLOID}" -eq 0 ]]; then
        echo "ERROR: No records remained after diploid filtering: \${DIPLOID_VCF}"
        exit 1
    fi

    echo "[phase] Checking for remaining non-diploid-style genotypes"
    N_BAD=\$(bcftools query -f '[%GT\n]' "\${DIPLOID_VCF}" | awk '
        \$1 !~ /^[.0-9]+[\/|][.0-9]+$/ {bad++}
        END {print bad+0}
    ')
    echo "[phase] Non-diploid-style GT records remaining: \${N_BAD}"

    if [[ "\${N_BAD}" -gt 0 ]]; then
        echo "ERROR: Non-diploid-style GT records still remain in \${DIPLOID_VCF}"
        exit 1
    fi

    echo "[phase] Running whatshap phase on diploid-only VCF"
    whatshap phase \
      --reference "\${REF}" \
      --ignore-read-groups \
      -o "\${PHASED_VCF}" \
      "\${DIPLOID_VCF}" "\${CLEAN_BAM}"

else
    echo "[phase] --external-phased-vcf supplied: \${EXTERNAL_PHASED_VCF}"
    echo "[phase] Skipping whatshap phase; using this externally-phased VCF instead."
    echo "[phase] (Phasing program is not restricted to whatshap -- haplotagging"
    echo "[phase]  and phase-block stats below work on any correctly phased VCF.)"

    if [[ ! -f "\${EXTERNAL_PHASED_VCF}" ]]; then
        echo "ERROR: external phased VCF not found: \${EXTERNAL_PHASED_VCF}"
        exit 1
    fi

    cp -f "\${EXTERNAL_PHASED_VCF}" "\${PHASED_VCF}"
fi

echo "[phase] Checking phased VCF is readable"
bcftools view -h "\${PHASED_VCF}" >/dev/null

echo "[phase] Rebuilding phased VCF index"
rm -f "\${PHASED_VCF}.csi" "\${PHASED_VCF}.tbi"
bcftools index -f "\${PHASED_VCF}"

echo "[phase] Verifying phased VCF carries FORMAT/PS tags for \${SAMPLE} (required by Part 2b onward)"
N_PS=\$(bcftools query -s "\${SAMPLE}" -f '[%PS\n]' "\${PHASED_VCF}" 2>/dev/null | awk '\$1!="." && \$1!=""' | wc -l)
echo "[phase] Records with a non-missing FORMAT/PS value: \${N_PS}"

if [[ "\${N_PS}" -eq 0 ]]; then
    echo "ERROR: \${PHASED_VCF} has no FORMAT/PS tags for \${SAMPLE}."
    echo "       Phasing did not produce phase sets, so Part 2b cannot proceed."
    echo "       If you supplied --external-phased-vcf, that tool's output must"
    echo "       include FORMAT/PS and phased (|) genotypes."
    exit 1
fi

echo "[phase] Writing phasing stats"
whatshap stats \
  --tsv="\${PHASE_STATS_TSV}" \
  --block-list="\${PHASE_BLOCKS_TSV}" \
  "\${PHASED_VCF}"

echo "[phase] Haplotagging cleaned BAM"
whatshap haplotag \
  --reference "\${REF}" \
  --ignore-read-groups \
  -o "\${HP_BAM}" \
  "\${PHASED_VCF}" \
  "\${CLEAN_BAM}"

echo "[phase] Indexing haplotagged BAM"
samtools index "\${HP_BAM}"

echo "[phase] Verifying haplotagged BAM integrity"
if ! samtools quickcheck -v "\${HP_BAM}"; then
    echo "ERROR: \${HP_BAM} fails quickcheck after haplotag"
    exit 1
fi

echo "[phase] Done"
EOF

    chmod +x "${OUT}"
    sbatch "${OUT}"
}


# ===========================================================================
## STEP 2b : Merge phased child VCF with unphased parent VCF, extract phased SNPs, and count shared alleles per haplotype block
## ==========================================================================


# ===========================================================================
## Helper: 
## Sometimes, after step 2b DNM window BED file is empty =>  0 LR-validated candidates
## So it returns 0. Will use in Parts 3, 3b, 4 to skip the local-rephasing refinement and just fix mutationn count = 0
## ==========================================================================
dnm_candidates_are_empty() {
    local MM_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis"
    local WIN_BED="${MM_DIR}/${SAMPLE_CHILD}_dnmc_plusminus${v}kb.bed"
    # Empty if the window BED is missing or has no data rows.
    if [[ ! -s "${WIN_BED}" ]]; then
        return 0
    fi
    return 1
}

# Still write a final_dnmc file that has the standard header but ZERO data rows, so
# Part 4's final_summary reports DNMC count = 0 instead of failing on a missing
# file. The header mirrors what count_mismatches.py emits for the dnmc.20kb.tsv.
write_empty_final_dnmc() {
    local PHASED_VCF="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local dnmc_file="${PHASED_VCF}/final_dnmc_${SAMPLE_CHILD}-from-${SAMPLE_PARENT}.tsv"
    if [[ ! -f "${dnmc_file}" ]]; then
        printf '# zero LR-validated DNM candidates\n' > "${dnmc_file}"
        echo "[zero-dnm] Wrote empty final DNM file: ${dnmc_file}"
    fi
}

# ===========================================================================
## Helper: ensure a bgzipped VCF has an index at least as new as the file.
## Rebuilds the .csi if it is missing OR older than the VCF, preventing the
## "[W::hts_idx_load3] The index file is older than the data file" warning and
## the inconsistent/partial reads it can cause. Safe to call repeatedly.
## ==========================================================================
ensure_fresh_vcf_index() {
    local vcf="$1"
    if [[ ! -f "${vcf}" ]]; then
        echo "[index] WARNING: VCF not found, cannot index: ${vcf}"
        return 0
    fi

    local csi="${vcf}.csi"
    local tbi="${vcf}.tbi"

    # Reindex if no index exists, or if either index is older than the VCF.
    if [[ ! -f "${csi}" && ! -f "${tbi}" ]] \
       || { [[ -f "${csi}" && "${vcf}" -nt "${csi}" ]]; } \
       || { [[ -f "${tbi}" && "${vcf}" -nt "${tbi}" ]]; }; then
        echo "[index] (Re)building index for ${vcf}"
        rm -f "${csi}" "${tbi}"
        bcftools index -f "${vcf}"
    else
        echo "[index] Index is up to date for ${vcf}"
    fi
}


# ===========================================================================
## Helper: gate Part 2b on Part 2's output actually being phased.
## Since users can have different phasing tools,
## It's better to check for the presence of FORMAT/PS and phased (|) genotypes in the
## phased child VCF before proceeding to Part 2b, rather than assuming whatshap
## Idea : Runs no matter HOW the phased child VCF got there -- whatshap phase (Part
## 2 default), --external-phased-vcf (Part 2, alternate phasing tool), or a
## phased VCF dropped into place by hand outside this script entirely. Part
## 2b/3/3b/4 all assume FORMAT/PS + phased (|) genotypes are present; without
## this gate a missing PS tag fails silently deep inside extract_phased_snp
## (zero rows extracted) instead of with a clear message right here.
## ==========================================================================
check_phased_vcf_has_ps() {
    local phased_vcf="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/${SAMPLE_CHILD}.illumVCF_LRbam.phased.vcf.gz"

    echo "[ps-check] Verifying Part 2 produced a phased VCF (FORMAT/PS present): ${phased_vcf}"

    if [[ ! -f "${phased_vcf}" ]]; then
        echo "ERROR: phased VCF not found: ${phased_vcf}"
        echo "       Run --part 2 first (with whatshap, or with --external-phased-vcf"
        echo "       pointing at output from another phasing program), or place a"
        echo "       phased VCF at this exact path yourself before running Part 2b."
        return 1
    fi

    if [[ ! -f "${phased_vcf}.csi" && ! -f "${phased_vcf}.tbi" ]]; then
        echo "[ps-check] No index found; indexing"
        bcftools index -f "${phased_vcf}"
    fi

    local n_ps
    n_ps=$(bcftools query -s "${SAMPLE_CHILD}" -f '[%PS\n]' "${phased_vcf}" 2>/dev/null \
        | awk '$1!="." && $1!=""' | wc -l)

    echo "[ps-check] Records with a non-missing FORMAT/PS value for ${SAMPLE_CHILD}: ${n_ps}"

    if [[ "${n_ps}" -eq 0 ]]; then
        echo "ERROR: ${phased_vcf} has no FORMAT/PS tags for ${SAMPLE_CHILD}."
        echo "       Part 2b requires a phased VCF with phase-set (PS) annotations."
        echo "       The phasing step itself is NOT restricted to whatshap -- any tool's"
        echo "       output is accepted as long as it writes FORMAT/PS and phased (|)"
        echo "       genotypes. Re-run Part 2 (optionally with --external-phased-vcf),"
        echo "       or regenerate the phased VCF with your tool of choice."
        return 1
    fi

    echo "[ps-check] OK: phased VCF carries PS tags for ${SAMPLE_CHILD}."
}



## ==========================================================================
## STEP 2b - Merge unphased parent VCF with phased child VCF,
## extract phased SNPs, and count shared alleles per haplotype block
## ==========================================================================
merge_unphased-parent_phased-child_vcfs() {
    # clean up
    rm -f ./slurm*out
    rm -f ./typescript

    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged

    mkdir -p ${PHASED_VCF}
    mkdir -p ${MERGED_PHASED_VCF}

    local PARENT_VCF=${ILLUM_DIR}/${SAMPLE_PARENT}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz
    local CHILD_VCF=${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_LRbam.phased.vcf.gz

    # Make sure indexes exist AND are newer than their VCFs. The old code only
    # indexed when an index was ABSENT, so a stale (older-than-VCF) index

    ensure_fresh_vcf_index "${PARENT_VCF}"
    ensure_fresh_vcf_index "${CHILD_VCF}"

    echo "merging samples"
    bcftools merge -m none -Oz \
        -o ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz \
        "${PARENT_VCF}" \
        "${CHILD_VCF}"

    local MERGED_VCF=${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz
    bcftools index -f "${MERGED_VCF}"
    echo "done"
}


## ==========================================================================
## Once the merged VCF is created, extract the phased SNPs and write to a TSV
## ==========================================================================
extract_phased_snp() {
    if command -v module &> /dev/null; then
        module load conda || true
    fi

    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged
    local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks

    mkdir -p ${EXTRACT_HAPLBLOCK}
    mkdir -p ${MERGED_PHASED_VCF}
    mkdir -p ${EXTRACT_HAPLBLOCK}

    ## Extract all the positions where we have a haplotype PS tag for the child
    bcftools view ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz \
    -i 'TYPE="snp" && N_ALT=1 && GT[0]!="mis" && GT[1]!="mis" && FMT/PS[1]!="."' \
    -Ou \
    | bcftools query \
        -s ${SAMPLE_PARENT},${SAMPLE_CHILD} \
        -f '%CHROM,%POS,[%PS,][%GT,][%DP,][%GQ,][%AD{0},][%AD{1},]\n' \
        > ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps.tsv


    ## Fix the extract file to tab delimited format
    conda activate "${CONDA_ENV_NAME}"

    python ${WORKING_DIR}/src/fix_PhaseSet.py \
      ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps.tsv \
      ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps_fixed.tsv
    
    # overwrite original file
    cp  ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps_fixed.tsv \
    ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps.tsv

    rm ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps_fixed.tsv 
    echo "Done fixing phase set file"
    echo "See file (below) for extracted snps"
    echo ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps.tsv
    }

## ==========================================================================
## Run the count_mismatches.py script to count shared alleles per haplotype block
## ==========================================================================
count_shared_alleles_per_PS_block() {
  local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
  local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks
  local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis

  mkdir -p ${MISMATCH_ANALYSIS}
  
  echo " Counting alleles per haplotype block"
  conda activate "${CONDA_ENV_NAME}"
  python ${WORKING_DIR}/src/count_mismatches.py \
    ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps.tsv \
    ${MISMATCH_ANALYSIS} \
    ${MIN_RDEPTH} \
    ${MAX_RDEPTH} \
    ${SAMPLE_PARENT} \
    ${SAMPLE_CHILD} \
    ${GT_QUAL} \
    ${NV_QUANTILE} \
    ${MM_DIFF_MIN} \
    0 ${NOTRECOUNT}
  echo "See outputs in ${MISMATCH_ANALYSIS}"
  echo "Done"
}


## ==========================================================================
## Filter the de novo mutation candidates based on genotype, depth, and allele-balance criteria
## ==========================================================================

filter_dnm_candidates() {
    echo "[DNM] Filtering de novo mutation candidates"

    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
    local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged

    local BED=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.bed
    local MERGED_VCF=${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz
    local OUT_VCF=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.vcf.gz
    local OUT_TSV=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.filt.tsv
    local OUT_BED=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.filt.bed
    local OUT_TSV_ORG=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.tsv
    local OUT_BED_ORG=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.bed

    if [[ ! -s "${BED}" ]]; then
        echo "ERROR: candidate BED not found or empty: ${BED}"
        return 1
    fi

    if [[ ! -f "${MERGED_VCF}" ]]; then
        echo "ERROR: merged VCF not found: ${MERGED_VCF}"
        return 1
    fi

    echo "[DNM] Subsetting merged VCF using BED regions (biallelic SNPs only)"

    # Restrict to biallelic SNPs in the candidate regions
    bcftools view \
        -R "${BED}" \
        -m2 -M2 -v snps \
        -Oz \
        -o "${OUT_VCF}" \
        "${MERGED_VCF}"

    bcftools index -f "${OUT_VCF}"

    # Sanity check: parent must be column 0, child column 1
    local col0
    col0=$(bcftools query -l "${OUT_VCF}" | sed -n '1p')
    if [[ "${col0}" != "${SAMPLE_PARENT}" ]]; then
        echo "ERROR: expected ${SAMPLE_PARENT} as sample[0] in ${MERGED_VCF}, got ${col0}."
        return 1
    fi

    echo "[DNM] Applying genotype / depth / allele-balance filters"

    # [0] = parent : must be hom-ref, no ALT-supporting reads
    # [1] = child  : must be phased het, with both REF and ALT support
    local filter_expr="
    (FMT/GT[1]==\"0|1\" || FMT/GT[1]==\"1|0\") &&
    FMT/GT[0]==\"0/0\" &&
    FMT/GT[0]!=\"1/1\" &&
    FMT/GT[0]!=\"0/1\" &&
    FMT/GT[0]!=\"1/0\" &&
    FMT/DP[0]>=${MIN_RDEPTH} && FMT/DP[0]<=${MAX_RDEPTH} &&
    FMT/DP[1]>=${MIN_RDEPTH} && FMT/DP[1]<=${MAX_RDEPTH} &&
    FMT/GQ[0]>=${GT_QUAL} &&
    FMT/GQ[1]>=${GT_QUAL} &&
    FMT/AD[1:1]>5 &&
    FMT/AD[1:0]>5 &&
    FMT/AD[0:1]==0
    "

    bcftools view "${OUT_VCF}" -i "${filter_expr}" \
    | bcftools query -H \
        -f '%CHROM\t%POS\t[%GT\t][%DP\t][%GQ\t][%AD\t]\n' \
    > "${OUT_TSV}"

    echo "[DNM] Creating filtered DNMC BED file"

    bcftools view "${OUT_VCF}" -i "${filter_expr}" \
    | bcftools query -f '%CHROM\t%POS\n' \
    | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' \
    > "${OUT_BED}"

    # Overwrite the unfiltered candidate files with the filtered set
    mv "${OUT_BED}" "${OUT_BED_ORG}"
    mv "${OUT_TSV}" "${OUT_TSV_ORG}"

    echo "[DNM] Done"
    echo "[DNM] Output VCF : ${OUT_VCF}"
    echo "[DNM] Output TSV : ${OUT_TSV_ORG}"
    echo "[DNM] Output BED : ${OUT_BED_ORG}"
    echo "[DNM] Candidate count: $(wc -l < "${OUT_BED_ORG}")"
}


## ==========================================================================
## Validate DNM candidates with long-read evidence
## ==========================================================================
validate_dnmc_with_long_reads(){

  conda activate "${CONDA_ENV_NAME}"

  local HIFI_BAM="${HP_BAM_PATH}"
  local OUT_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis"
  local DNM_BED="${OUT_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.bed"
  echo "[dnmc-check] Validating DNM candidates with long-read evidence"
  echo "[dnmc-check] LR BAM: ${HIFI_BAM}"

  python "${WORKING_DIR}/src/dnmc_readcheck.py" \
    "${SAMPLE_CHILD}" \
    "${HIFI_BAM}" \
    "${DNM_BED}" \
    "${MIN_BASE_QUAL}" \
    "${MIN_MAP_QUAL}" \
    "${WINDW}" \
    "${ALT_READ_COUNT}" \
    "T" \
    "dnmc" \
    "${TOTAL_READ_COUNT_MIN}" \
    "${OUT_DIR}"

    echo "[dnmc-check] Done validating DNMC candidates with long-read evidence"

}


## ==========================================================================
##  Additional validation of DNM candidates with parent Illumina BAM
## ==========================================================================
validate_dnmc_with_parent_bam() {
    echo "[parent-check] Validating DNM candidates against parent Illumina BAM"

    conda activate "${CONDA_ENV_NAME}"

    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
    local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged

    local IN_BED="${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.bed"
    local OUT_BED="${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.parentBAM.bed"
    local MERGED_VCF="${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz"
    local PARENT_BAM="${ORIG_BAM_PARENT_PATH}"
    local MAX_PARENT_ALT=1

    if [[ ! -s "${IN_BED}" ]]; then
        echo "[parent-check] No candidates in ${IN_BED}; nothing to do"
        return 0
    fi

    if [[ ! -f "${PARENT_BAM}" ]]; then
        echo "[parent-check] WARNING: parent BAM not found at ${PARENT_BAM}"
        echo "[parent-check] Skipping parent-BAM validation (results may include FN-parent DNMs)"
        return 0
    fi

    # Make sure the parent BAM is indexed
    if [[ ! -f "${PARENT_BAM}.bai" && ! -f "${PARENT_BAM%.bam}.bai" ]]; then
        echo "[parent-check] Indexing parent BAM"
        samtools index "${PARENT_BAM}"
    fi

    local N_BEFORE
    N_BEFORE=$(wc -l < "${IN_BED}")

    python "${WORKING_DIR}/src/parent_readcheck.py" \
        "${SAMPLE_PARENT}" \
        "${PARENT_BAM}" \
        "${MERGED_VCF}" \
        "${IN_BED}" \
        "${OUT_BED}" \
        "${MIN_BASE_QUAL}" \
        "${MIN_MAP_QUAL}" \
        "${MAX_PARENT_ALT}" \
        "${MISMATCH_ANALYSIS}"

    # Swap filtered BED into the canonical location
    mv "${OUT_BED}" "${IN_BED}"

    local N_AFTER
    N_AFTER=$(wc -l < "${IN_BED}")

    echo "[parent-check] Candidates before: ${N_BEFORE}"
    echo "[parent-check] Candidates after : ${N_AFTER}"
    echo "[parent-check] Removed by parent-BAM evidence: $((N_BEFORE - N_AFTER))"
}




# ===========================================================================
## STEP 3: Local rephasing around candidate DNMs  (Option-A refactor)
## ==========================================================================

_run_local_phasing_body() {
    set -euo pipefail

    if command -v module &> /dev/null; then
        module load bcftools || true
        module load samtools || true
        module load conda || true
    fi
    conda activate "${CONDA_ENV_NAME}"

    local SAMPLE="${SAMPLE_CHILD}"
    local MM_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis/"

    local ILLUM_VCF="${ILLUM_DIR}/${SAMPLE}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz"
    local CLEAN_BAM="${CLEAN_BAM_PATH}"

    local REGIONS_BED="${MM_DIR}/${SAMPLE}_dnmc_plusminus${v}kb.bed"

    mkdir -p "${MM_DIR}"
    cd "${MM_DIR}"

    local LOCAL_VCF="${SAMPLE}_dnmc_plusminus${v}kb.vcf.gz"
    local LOCAL_DIPLOID_VCF="${SAMPLE}_dnmc_plusminus${v}kb.diploid.vcf.gz"
    local LOCAL_PHASED_VCF="${SAMPLE}.illumVCF_LRbam.phased.${v}kb.vcf.gz"
    local LOCAL_STATS_TSV="${SAMPLE}.illumVCF_LRbam.stats.${v}kb.tsv"
    local LOCAL_BLOCKS_TSV="${SAMPLE}.illumVCF_LRbam.blocks.${v}kb.tsv"

    echo "[local phase] Regions BED        : ${REGIONS_BED}"
    echo "[local phase] Working dir        : ${MM_DIR}"

    if [[ ! -f "${ILLUM_VCF}" ]]; then echo "ERROR: Illumina VCF not found: ${ILLUM_VCF}"; return 1; fi
    if [[ ! -f "${CLEAN_BAM}" ]]; then echo "ERROR: Clean BAM not found: ${CLEAN_BAM}"; return 1; fi
    if [[ ! -f "${REGIONS_BED}" ]]; then echo "ERROR: Regions BED not found: ${REGIONS_BED}"; return 1; fi
    if [[ ! -f "${CLEAN_BAM}.bai" ]]; then samtools index "${CLEAN_BAM}"; fi

    rm -f "${LOCAL_VCF}" "${LOCAL_VCF}.csi" "${LOCAL_VCF}.tbi"
    rm -f "${LOCAL_DIPLOID_VCF}" "${LOCAL_DIPLOID_VCF}.csi" "${LOCAL_DIPLOID_VCF}.tbi"
    rm -f "${LOCAL_PHASED_VCF}" "${LOCAL_PHASED_VCF}.csi" "${LOCAL_PHASED_VCF}.tbi"
    rm -f "${LOCAL_STATS_TSV}" "${LOCAL_BLOCKS_TSV}"

    echo "[local phase] Subsetting Illumina VCF to local regions"
    bcftools view -R "${REGIONS_BED}" -s "${SAMPLE}" -Oz -o "${LOCAL_VCF}" "${ILLUM_VCF}"
    bcftools index -f "${LOCAL_VCF}"

    local N_LOCAL
    N_LOCAL=$(bcftools index -n "${LOCAL_VCF}")
    if [[ "${N_LOCAL}" -eq 0 ]]; then echo "ERROR: No variants found in local regions"; return 1; fi

    echo "[local phase] Building diploid-only local VCF"
    bcftools view -Ou "${LOCAL_VCF}" \
      | bcftools view -i 'GT="het" || GT="hom" || GT="ref" || GT="alt"' \
      -Oz -o "${LOCAL_DIPLOID_VCF}"
    bcftools index -f "${LOCAL_DIPLOID_VCF}"

    local N_DIPLOID
    N_DIPLOID=$(bcftools index -n "${LOCAL_DIPLOID_VCF}")
    if [[ "${N_DIPLOID}" -eq 0 ]]; then echo "ERROR: No diploid variants remained"; return 1; fi

    echo "[local phase] Running whatshap"
    whatshap phase --reference "${REF}" --ignore-read-groups \
      -o "${LOCAL_PHASED_VCF}" "${LOCAL_DIPLOID_VCF}" "${CLEAN_BAM}"

    bcftools view -h "${LOCAL_PHASED_VCF}" >/dev/null
    rm -f "${LOCAL_PHASED_VCF}.csi" "${LOCAL_PHASED_VCF}.tbi"
    bcftools index -f "${LOCAL_PHASED_VCF}"

    whatshap stats --tsv="${LOCAL_STATS_TSV}" --block-list="${LOCAL_BLOCKS_TSV}" "${LOCAL_PHASED_VCF}"
    echo "[local phase] Done (inline)."
}

regenerate_phasing_job() {
    local OUT="${PRJ_DIR}/build_hapl_${SAMPLE_CHILD}.slurm"
    local SCRIPT_PATH
    SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
    local COMMON_ARGS
    COMMON_ARGS="$(_chain_common_args)"

    cat << SLURM_EOF > "${OUT}"
#!/bin/bash
#SBATCH -J build_hapl_${SAMPLE_CHILD}_localphase
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME}
#SBATCH --mem=8G
#SBATCH -A ${ACCOUNT}
${SBATCH_PARTITION_LINE}
#SBATCH -o slurm_output/build_hapl_local_%j.txt
#SBATCH -e slurm_output/build_hapl_local_%j.err

set -euo pipefail
echo "[local phase job] Start: \$(date)"
bash "${SCRIPT_PATH}" ${COMMON_ARGS} --part 3-inline
echo "[local phase job] Done: \$(date)"
SLURM_EOF

    chmod +x "${OUT}"
    mkdir -p slurm_output
    sbatch "${OUT}"
}

# ===========================================================================
## CHAIN : Submit Parts 2b -> 3 -> 3b -> 4 as four Slurm jobs chained with
##         afterok dependencies. 
## ==========================================================================

_chain_common_args() {
    printf '%s' \
"--prj-dir ${PRJ_DIR} \
--sample-child ${SAMPLE_CHILD} \
--sample-parent ${SAMPLE_PARENT} \
--source $(basename "${SOURCE_FILE}") \
--reference ${REF} \
--cpus ${CPUS} \
--time ${TIME} \
--min-rdepth ${MIN_RDEPTH} \
--max-rdepth ${MAX_RDEPTH} \
--gt-qual ${GT_QUAL} \
--nv-quantile ${NV_QUANTILE} \
--mm-diff-min ${MM_DIFF_MIN} \
--min-base-qual ${MIN_BASE_QUAL} \
--min-map-qual ${MIN_MAP_QUAL} \
--window ${WINDW} \
--alt-read-count ${ALT_READ_COUNT} \
--verbose ${VERBOSE} \
--exclude-chroms ${EXCLUDE_CHROMS} \
--total-rd-ct-min ${TOTAL_READ_COUNT_MIN} \
--window_rephase ${WINDOW_REPHASE} \
--threshold_rephase ${THRESHOLD_REPHASE} \
--conda-env ${CONDA_ENV_NAME}"
}

_chain_write_step_slurm() {
    local part_token="$1" out_slurm="$2" jobname="$3" walltime="$4" mem="$5"
    local log_dir="$6" common_args="$7" script_path="$8"

    cat << STEP_EOF > "${out_slurm}"
#!/bin/bash
#SBATCH -J ${jobname}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${walltime}
#SBATCH --mem=${mem}
#SBATCH -A ${ACCOUNT}
${SBATCH_PARTITION_LINE}
#SBATCH -o ${log_dir}/${jobname}_%j.out
#SBATCH -e ${log_dir}/${jobname}_%j.err

set -euo pipefail
echo "[chain:${part_token}] Start: \$(date)"
bash "${script_path}" ${common_args} --part ${part_token}
echo "[chain:${part_token}] Done: \$(date)"
STEP_EOF
    chmod +x "${out_slurm}"
}

submit_chain_2b_4_jobs() {
    echo "[chain] Submitting 2b -> 3 -> 3b -> 4 as an afterok Slurm chain"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local CHAIN_DIR="${PHASED_DIR}/chain_2b_4"
    local LOG_DIR="${CHAIN_DIR}/logs"
    local SCRIPT_PATH
    SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

    mkdir -p "${CHAIN_DIR}" "${LOG_DIR}" slurm_output

    local COMMON_ARGS
    COMMON_ARGS="$(_chain_common_args)"

    local S2B="${CHAIN_DIR}/run_2b_${SAMPLE_CHILD}.slurm"
    _chain_write_step_slurm "2b" "${S2B}" "chain2b_${SAMPLE_CHILD}" "${TIME}" "12G" "${LOG_DIR}" "${COMMON_ARGS}" "${SCRIPT_PATH}"
    local J2B; J2B=$(sbatch --parsable "${S2B}")
    echo "[chain] Submitted 2b : ${J2B}"

    local S3="${CHAIN_DIR}/run_3_${SAMPLE_CHILD}.slurm"
    _chain_write_step_slurm "3-inline" "${S3}" "chain3_${SAMPLE_CHILD}" "${TIME}" "12G" "${LOG_DIR}" "${COMMON_ARGS}" "${SCRIPT_PATH}"
    local J3; J3=$(sbatch --parsable --dependency=afterok:${J2B} "${S3}")
    echo "[chain] Submitted 3  : ${J3} (after ${J2B})"

    local S3B="${CHAIN_DIR}/run_3b_${SAMPLE_CHILD}.slurm"
    _chain_write_step_slurm "3b" "${S3B}" "chain3b_${SAMPLE_CHILD}" "${TIME}" "12G" "${LOG_DIR}" "${COMMON_ARGS}" "${SCRIPT_PATH}"
    local J3B; J3B=$(sbatch --parsable --dependency=afterok:${J3} "${S3B}")
    echo "[chain] Submitted 3b : ${J3B} (after ${J3})"

    local S4="${CHAIN_DIR}/run_4_${SAMPLE_CHILD}.slurm"
    _chain_write_step_slurm "4" "${S4}" "chain4_${SAMPLE_CHILD}" "${TIME}" "12G" "${LOG_DIR}" "${COMMON_ARGS}" "${SCRIPT_PATH}"
    local J4; J4=$(sbatch --parsable --dependency=afterok:${J3B} "${S4}")
    echo "[chain] Submitted 4  : ${J4} (after ${J3B})"

    echo "[chain] Chain: 2b(${J2B}) -> 3(${J3}) -> 3b(${J3B}) -> 4(${J4})"
    echo "[chain] Logs in: ${LOG_DIR}"
}

# ===========================================================================
## STEP 3b / PART 3b : Re-merge phased child VCF with unphased parent VCF, extract phased SNPs, and recount shared alleles per haplotype block based on new phasing
## ==========================================================================

remerge_unphased-parent_phased-child_vcfs() {
    echo "Re-merging parent vcf and child vcf"
    rm -f ./slurm*out
    rm -f ./typescript

    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/
    local MM_DIR=${PHASED_VCF}/mismatch_analysis
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged

    mkdir -p ${PHASED_VCF}
    mkdir -p ${MERGED_PHASED_VCF}

    # move child vcf to correct place
    ls ${MM_DIR}/${SAMPLE_CHILD}.illumVCF_LRbam* 1>/dev/null 2>&1 && \
        mv -f ${MM_DIR}/${SAMPLE_CHILD}.illumVCF_LRbam* ${PHASED_VCF}/

    local PARENT_VCF=${ILLUM_DIR}/${SAMPLE_PARENT}.CHM13.illumina.unphased.noPS.vcf.gz
    local CHILD_VCF=${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_LRbam.phased.${v}kb.vcf.gz

    # Freshness-aware indexing (rebuilds if missing OR older than the VCF).
    ensure_fresh_vcf_index "${PARENT_VCF}"
    ensure_fresh_vcf_index "${CHILD_VCF}"

    bcftools merge -m none -Oz \
        -o ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.${v}kb.vcf.gz \
        "${PARENT_VCF}" \
        "${CHILD_VCF}"

    ensure_fresh_vcf_index "${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.${v}kb.vcf.gz"
    echo "done"
}

reextract_phased_snp() {
    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged
    local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks

    if command -v module &> /dev/null; then
        module load bcftools || true
        module load conda || true
    fi

    echo "Reextracting phased snps"
    

    mkdir -p ${EXTRACT_HAPLBLOCK}
    mkdir -p ${MERGED_PHASED_VCF}
    mkdir -p ${EXTRACT_HAPLBLOCK}

    bcftools view ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.${v}kb.vcf.gz \
    -i 'TYPE="snp" && N_ALT=1 && GT[0]!="mis" && GT[1]!="mis" && FMT/PS[1]!="."' \
    -Ou \
    | bcftools query \
        -s ${SAMPLE_PARENT},${SAMPLE_CHILD} \
        -f '%CHROM,%POS,[%PS,][%GT,][%DP,][%GQ,][%AD{0},][%AD{1},]\n' \
        > ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps.tsv

 
    ## fix again the extracted format into tab delim file
    conda activate "${CONDA_ENV_NAME}"
    python ${WORKING_DIR}/src/fix_PhaseSet.py \
      ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps.tsv \
      ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps_fixed.tsv
    
    # overwrite original file
    cp  ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps_fixed.tsv \
    ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps.tsv

    rm ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps_fixed.tsv 
    echo "Done fixing phase set file"
    echo "See file (below) for extracted snps"
    echo ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps.tsv
    }



recount_shared_alleles_per_PS_block() {
  local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
  local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks
  local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis
  local RECOUNT="T"
  local DNMC_FILE=${PHASED_VCF}/mismatch_analysis/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.tsv

  echo "Recounting mismatches in PS blocks"
  
  conda activate "${CONDA_ENV_NAME}"
  python "${WORKING_DIR}/src/count_mismatches.py" \
  "${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps.tsv" \
  "${MISMATCH_ANALYSIS}" \
  "${MIN_RDEPTH}" \
  "${MAX_RDEPTH}" \
  "${SAMPLE_PARENT}" \
  "${SAMPLE_CHILD}" \
  "${GT_QUAL}" \
    0 \
    "${MM_DIFF_MIN}" \
    "${WINDW}" \
    "${RECOUNT}" \
    "${DNMC_FILE}"

    echo "See outputs in ${MISMATCH_ANALYSIS}"
    echo "Done"

}

clean_up(){
  local MM_DIR=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis
  local SLICE_DIR=${MM_DIR}/sliced/

  mkdir -p ${SLICE_DIR}

  # After the output-path migration, Part 2b / dnmc_readcheck.py writes the
  # *_plusminus*.bed and *_LR_validated_* files into mismatch_analysis/, not the
  # current working directory. Move whatever is present from there, and never
  # fail the job if a class of file is absent (nullglob + guarded mv).
  shopt -s nullglob
  local f
  for f in "${MM_DIR}"/${SAMPLE_CHILD}*plusminus* "${MM_DIR}"/${SAMPLE_CHILD}*LR_validated*; do
      mv -f "${f}" "${SLICE_DIR}/"
  done
  # Also sweep any stragglers still left in CWD by older/manual runs (guarded).
  for f in ./${SAMPLE_CHILD}*plusminus* ./${SAMPLE_CHILD}*LR_validated*; do
      mv -f "${f}" "${SLICE_DIR}/"
  done
  shopt -u nullglob
}


calculate_callable_genome() {
  local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
  local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks
  local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis
  local POST_ANALYSIS=${PHASED_VCF}/mismatch_analysis/denum_calcul
  
  mkdir -p ${POST_ANALYSIS}
  conda activate "${CONDA_ENV_NAME}"
  EXCLUDE_CHROMS="${EXCLUDE_CHROMS}" \
  
  python "${WORKING_DIR}/src/callable_genome.py" \
  "${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps.tsv" \
  "${POST_ANALYSIS}" \
  "${MIN_RDEPTH}" \
  "${MAX_RDEPTH}" \
  "${SAMPLE_PARENT}" \
  "${SAMPLE_CHILD}" \
  "${GT_QUAL}" \
  "${NV_QUANTILE}" \
  "${MM_DIFF_MIN}" \
  0 "${NOTRECOUNT}"
  
  echo "See outputs in ${MISMATCH_ANALYSIS}"
  echo "Done"

}



REfilter_dnm_candidates() {
    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
    local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis/denum_calcul
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged
    local BED=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.bed
    local MERGED_VCF=${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz
    local OUT_VCF=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.vcf.gz
    local OUT_TSV=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.filt.tsv
    local OUT_BED=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.filt.bed
    local OUT_TSV_ORG=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.tsv
    local OUT_BED_ORG=${MISMATCH_ANALYSIS}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.bed


    echo "[DNM] Subsetting merged VCF using BED regions (biallelic SNPs only)"

    # Restrict to biallelic SNPs
    bcftools view \
        -R "${BED}" \
        -m2 -M2 -v snps \
        -Oz \
        -o "${OUT_VCF}" \
        "${MERGED_VCF}"

    bcftools index -f "${OUT_VCF}"

    echo "[DNM] Applying genotype / depth / allele-balance filters"

    filter_expr="
    (FMT/GT[1]==\"0|1\" || FMT/GT[1]==\"1|0\") &&
    (FMT/GT[0]==\"0/1\" || FMT/GT[0]==\"1/0\") &&

    FMT/DP>=${MIN_RDEPTH} && FMT/DP<=${MAX_RDEPTH} &&
    FMT/GQ>=${GT_QUAL} &&

    FMT/AD[1:1]>5 &&
    FMT/AD[1:0]>5
    "

    bcftools view "${OUT_VCF}" -i "${filter_expr}" \
    | bcftools query -H \
        -f '%CHROM\t%POS\t[%GT\t][%DP\t][%GQ\t][%AD\t]\n' \
    > "${OUT_TSV}"

    # Create a bed file
    echo "[DNM] Creating filtered DNMC BED file"

    bcftools view "${OUT_VCF}" -i "${filter_expr}" \
    | bcftools query -f '%CHROM\t%POS\n' \
    | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' \
    > "${OUT_BED}"

    # Clean up
    mv ${OUT_BED} ${OUT_BED_ORG}
    mv ${OUT_TSV} ${OUT_TSV_ORG}

    echo "[DNM] Done"
    echo "[DNM] Output VCF: ${OUT_VCF}"
    echo "[DNM] Output TSV: ${OUT_TSV_ORG}"
    echo "[DNM] Output BED: ${OUT_BED_ORG}"

}


REvalidate_dnmc_with_long_reads(){

  conda activate "${CONDA_ENV_NAME}"

  local OUT_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis/denum_calcul"
  local HIFI_BAM="${HP_BAM_PATH}"
  local DNM_BED="${OUT_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.bed"
  local VERBOSE="F"

  python "${WORKING_DIR}/src/dnmc_readcheck.py" \
    "${SAMPLE_CHILD}" \
    "${HIFI_BAM}" \
    "${DNM_BED}" \
    "${MIN_BASE_QUAL}" \
    "${MIN_MAP_QUAL}" \
    "${WINDW}" \
    "${ALT_READ_COUNT}" \
    "${VERBOSE}" \
    "hetc" \
    "${TOTAL_READ_COUNT_MIN}" \
    "${OUT_DIR}"
}

# ===========================================================================
## STEP 4 : Summarize results, calculate final mutation rate, and generate final report
## ==========================================================================
final_summary() {

    local WORKING_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local DENUM_DIR="${WORKING_DIR}/mismatch_analysis/denum_calcul"
    local dnmc_file="${WORKING_DIR}/final_dnmc_${SAMPLE_CHILD}-from-${SAMPLE_PARENT}.tsv"
    local nb_qualified_snps="${DENUM_DIR}/${SAMPLE_CHILD}_LR_validated_hetc.bed"
    local callable_genome_file="${DENUM_DIR}/callable_genome.txt"

    mkdir -p "${DENUM_DIR}"

    # NOTE: dnmc_readcheck.py now writes its outputs directly into
    # ${DENUM_DIR} (the hetc *_LR_validated_*/*_plusminus* files) and into
    # mismatch_analysis (the dnmc files), so the old "compgen ... mv from CWD"
    # cleanup is no longer needed and has been removed. Only the dnmc.20kb.tsv
    # produced by Part 3b (count_mismatches.py) still needs relocating below.

    if [[ -f "${WORKING_DIR}/mismatch_analysis/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.20kb.tsv" ]]; then
        mv -f "${WORKING_DIR}/mismatch_analysis/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.20kb.tsv" \
              "${dnmc_file}"
    fi

    local dnmc_count
    dnmc_count=$(($(wc -l < "${dnmc_file}") - 1))

    local qualified_count
    qualified_count=$(wc -l < "${nb_qualified_snps}")

    local total_sampled
    total_sampled=$(grep "total_sampled_snps" "${callable_genome_file}" | awk '{print $2}')

    local accessible_bases
    accessible_bases=$(grep "total_callable_bases" "${callable_genome_file}" | awk '{print $2}')

    local mutation_rate="NA"

    if [[ -n "${total_sampled}" && -n "${accessible_bases}" && "${total_sampled}" -gt 0 ]]; then
        mutation_rate=$(awk -v d="${dnmc_count}" \
                            -v q="${qualified_count}" \
                            -v s="${total_sampled}" \
                            -v c="${accessible_bases}" \
                            'BEGIN { printf "%.6e", d / ((q/s) * c) }')
        callable_bases=$(awk -v q="${qualified_count}" \
                            -v s="${total_sampled}" \
                            -v c="${accessible_bases}" \
                            'BEGIN { printf "%.6e", (q/s) * c }')

    fi

    echo "================ Final results summary ================="
    echo "De novo candidate        : ${dnmc_file}"
    echo "Child sample             : ${SAMPLE_CHILD}"
    echo "Parent sample            : ${SAMPLE_PARENT}"
    echo "Min read depth (SR)      : ${MIN_RDEPTH}"
    echo "Max read depth (SR)      : ${MAX_RDEPTH}"
    echo "Genotype qual cutoff (SR): ${GT_QUAL}"
    echo "Noise quantile           : ${NV_QUANTILE}"
    echo "Mismatch diff threshold  : ${MM_DIFF_MIN}"
    echo "Min base quality (SR)    : ${MIN_BASE_QUAL}"
    echo "Min mapping quality (SR) : ${MIN_MAP_QUAL}"
    echo "Window size (bp)         : ${WINDW}"
    echo "Alt read count (LR)      : ${ALT_READ_COUNT}"
    echo "------------------------------------------------"
    echo "DNMC count               : ${dnmc_count}"
    echo "Number of snps qualified : ${qualified_count}"
    echo "Number of snps sampled   : ${total_sampled}"
    echo "Accessible genome size   : ${accessible_bases}"
    echo "Callable genome size     : ${callable_bases} [=(snp qualified / snp sampled) * accessible genome]"
    echo "Mutation rate            : ${mutation_rate}"
    echo "========================================================="
    echo
}



# ===========================================================================
## STEP 5 : Rephase blocks with many mismatches and re-evaluate candidate DNMs
## ==========================================================================

rephase_blocks_with_mismatches(){
    if command -v module &> /dev/null; then
        module load conda || true
    fi
    conda activate "${CONDA_ENV_NAME}"

    echo "[Rephasing]: Rephasing blocks with many mismatches (potential switch errors)"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local REPHASE_DIR="${PHASED_DIR}/rephased_blocks"
    local LOCAL_VCF_DIR="${REPHASE_DIR}/local_vcfs"
    local LOG_DIR="${REPHASE_DIR}/logs"

    local REPHASE_BED="${REPHASE_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_rephase_regions.bed"
    local FULL_VCF="${ILLUM_DIR}/${SAMPLE_CHILD}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz"
    local FULL_BAM="${CLEAN_BAM_PATH}"

    mkdir -p "${REPHASE_DIR}"
    mkdir -p "${LOCAL_VCF_DIR}"
    mkdir -p "${LOG_DIR}"

    rm -f "${REPHASE_BED}"

    echo "[Rephasing]: Detecting regions to rephase"

    python "${WORKING_DIR}/src/rephase_blocks.py" \
        --in_block_summary_file "${PHASED_DIR}/mismatch_analysis/${SAMPLE_PARENT}_${SAMPLE_CHILD}_mismatch.tsv" \
        --in_ht_file "${PHASED_DIR}/HTblocks/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps.tsv" \
        --out_dir "${REPHASE_DIR}" \
        --parent_id "${SAMPLE_PARENT}" \
        --child_id "${SAMPLE_CHILD}" \
        --window_size "${WINDOW_REPHASE}" \
        --threshold "${THRESHOLD_REPHASE}"

    if [[ ! -f "${REPHASE_BED}" ]]; then
        echo "[Rephasing]: No BED file produced by Python. Skipping local rephasing."
        return 0
    fi

    if [[ ! -s "${REPHASE_BED}" ]]; then
        echo "[Rephasing]: BED file exists but contains no regions. Skipping local rephasing."
        return 0
    fi

    echo "[Rephasing]: Regions written to ${REPHASE_BED}"
    echo "[Rephasing]: Keeping local phased VCFs in ${LOCAL_VCF_DIR}"

    if [[ ! -f "${FULL_VCF}" ]]; then
        echo "ERROR: child VCF not found: ${FULL_VCF}"
        return 1
    fi

    if [[ ! -f "${FULL_BAM}" ]]; then
        echo "ERROR: clean BAM not found: ${FULL_BAM}"
        return 1
    fi

    if [[ ! -f "${FULL_BAM}.bai" ]]; then
        echo "[Rephasing]: BAM index not found. Creating index."
        samtools index "${FULL_BAM}"
    fi

    local TOTAL_REGIONS
    TOTAL_REGIONS=$(grep -cv '^[[:space:]]*$' "${REPHASE_BED}")

    if [[ -z "${TOTAL_REGIONS}" || "${TOTAL_REGIONS}" -lt 1 ]]; then
        echo "[Rephasing]: No valid regions found in ${REPHASE_BED}. Skipping local rephasing."
        return 0
    fi

    echo "[Rephasing]: Total regions to submit: ${TOTAL_REGIONS}"

    local MAX_REPHASE_JOBS=200

    local N_ARRAY_TASKS
    if [[ "${TOTAL_REGIONS}" -lt "${MAX_REPHASE_JOBS}" ]]; then
        N_ARRAY_TASKS="${TOTAL_REGIONS}"
    else
        N_ARRAY_TASKS="${MAX_REPHASE_JOBS}"
    fi

    local ARRAY_SPEC="1-${N_ARRAY_TASKS}"

    echo "[Rephasing]: Total regions: ${TOTAL_REGIONS}"
    echo "[Rephasing]: Number of Slurm array tasks: ${N_ARRAY_TASKS}"
    echo "[Rephasing]: Array spec: ${ARRAY_SPEC}"

    local REPHASE_CPUS=1
    local REPHASE_MEM="8G"
    local REPHASE_TIME="04:00:00"

    local CONFIG_FILE="${REPHASE_DIR}/rephase_config.env"
    local ARRAY_SCRIPT="${REPHASE_DIR}/run_rephase_region_array.sh"
    local SUMMARY_SCRIPT="${REPHASE_DIR}/summarize_rephase_array.sh"

    echo "[Rephasing]: Removing old success/failure markers"
    rm -f "${LOCAL_VCF_DIR}"/region_*.success
    rm -f "${LOCAL_VCF_DIR}"/region_*.failed

    echo "[Rephasing]: Writing config to ${CONFIG_FILE}"

    cat > "${CONFIG_FILE}" <<EOF
SAMPLE_PARENT="${SAMPLE_PARENT}"
SAMPLE_CHILD="${SAMPLE_CHILD}"
REPHASE_BED="${REPHASE_BED}"
LOCAL_VCF_DIR="${LOCAL_VCF_DIR}"
FULL_VCF="${FULL_VCF}"
FULL_BAM="${FULL_BAM}"
REF="${REF}"
REPHASE_DIR="${REPHASE_DIR}"
CONDA_ENV_NAME="${CONDA_ENV_NAME}"
EOF

    echo "[Rephasing]: Writing Slurm array script to ${ARRAY_SCRIPT}"

    cat > "${ARRAY_SCRIPT}" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="$1"
N_ARRAY_TASKS="$2"

source "${CONFIG_FILE}"

if command -v module &> /dev/null; then
    module load conda || true
    module load samtools || true
    module load bcftools || true
fi

conda activate "${CONDA_ENV_NAME}"

TASK_ID="${SLURM_ARRAY_TASK_ID}"

echo "[Rephasing-array]: Task ${TASK_ID} of ${N_ARRAY_TASKS}"
echo "[Rephasing-array]: BED file: ${REPHASE_BED}"

awk -v task_id="${TASK_ID}" -v n_tasks="${N_ARRAY_TASKS}" '
    NF >= 3 {
        row += 1
        if (((row - 1) % n_tasks) == (task_id - 1)) {
            print row "\t" $0
        }
    }
' "${REPHASE_BED}" | while IFS=$'\t' read -r REGION_INDEX chrom start end rest; do

    REGION="${chrom}:${start}-${end}"

    echo "============================================================"
    echo "[Rephasing-array]: Task ${TASK_ID}, BED row ${REGION_INDEX}"
    echo "[Rephasing-array]: Region ${REGION}"

    REGION_PREFIX="${LOCAL_VCF_DIR}/region_${REGION_INDEX}_${chrom}_${start}_${end}"
    REGION_VCF="${REGION_PREFIX}.vcf.gz"
    REGION_BAM="${REGION_PREFIX}.bam"
    REGION_PHASED_VCF="${REGION_PREFIX}.phased.vcf.gz"
    REGION_HAPLOTAGGED_BAM="${REGION_PREFIX}.haplotagged.bam"
    SUCCESS_MARKER="${REGION_PREFIX}.success"
    FAIL_MARKER="${REGION_PREFIX}.failed"

    rm -f "${SUCCESS_MARKER}" "${FAIL_MARKER}"

    rm -f "${REGION_VCF}" "${REGION_VCF}.tbi" "${REGION_VCF}.csi"
    rm -f "${REGION_BAM}" "${REGION_BAM}.bai"
    rm -f "${REGION_PHASED_VCF}" "${REGION_PHASED_VCF}.tbi" "${REGION_PHASED_VCF}.csi"
    rm -f "${REGION_HAPLOTAGGED_BAM}" "${REGION_HAPLOTAGGED_BAM}.bai"

    echo "[Rephasing-array]: Extracting region from original child VCF"

    if ! bcftools view \
        -r "${REGION}" \
        -s "${SAMPLE_CHILD}" \
        -Oz \
        -o "${REGION_VCF}" \
        "${FULL_VCF}"
    then
        echo "[Rephasing-array]: FAILED bcftools view for ${REGION}"
        touch "${FAIL_MARKER}"
        continue
    fi

    if ! bcftools index -f "${REGION_VCF}"; then
        echo "[Rephasing-array]: FAILED bcftools index for ${REGION_VCF}"
        touch "${FAIL_MARKER}"
        continue
    fi

    NVAR=$(bcftools view -H "${REGION_VCF}" | wc -l)

    if [[ "${NVAR}" -eq 0 ]]; then
        echo "[Rephasing-array]: No variants in ${REGION}; skipping."
        touch "${SUCCESS_MARKER}"
        continue
    fi

    echo "[Rephasing-array]: Extracting reads overlapping region"

    if ! samtools view \
        -b \
        "${FULL_BAM}" \
        "${REGION}" \
        -o "${REGION_BAM}"
    then
        echo "[Rephasing-array]: FAILED samtools view for ${REGION}"
        touch "${FAIL_MARKER}"
        continue
    fi

    if ! samtools index "${REGION_BAM}"; then
        echo "[Rephasing-array]: FAILED samtools index for ${REGION_BAM}"
        touch "${FAIL_MARKER}"
        continue
    fi

    echo "[Rephasing-array]: Running whatshap phase"

    if ! whatshap phase \
        --reference "${REF}" \
        --ignore-read-groups \
        --sample "${SAMPLE_CHILD}" \
        -o "${REGION_PHASED_VCF}" \
        "${REGION_VCF}" \
        "${REGION_BAM}"
    then
        echo "[Rephasing-array]: FAILED whatshap phase for ${REGION}"
        touch "${FAIL_MARKER}"
        continue
    fi

    if ! bcftools index -f "${REGION_PHASED_VCF}"; then
        echo "[Rephasing-array]: FAILED bcftools index for ${REGION_PHASED_VCF}"
        touch "${FAIL_MARKER}"
        continue
    fi

    echo "[Rephasing-array]: Running whatshap haplotag"

    if ! whatshap haplotag \
        --reference "${REF}" \
        --ignore-read-groups \
        --sample "${SAMPLE_CHILD}" \
        -o "${REGION_HAPLOTAGGED_BAM}" \
        "${REGION_PHASED_VCF}" \
        "${REGION_BAM}"
    then
        echo "[Rephasing-array]: FAILED whatshap haplotag for ${REGION}"
        touch "${FAIL_MARKER}"
        continue
    fi

    if ! samtools index "${REGION_HAPLOTAGGED_BAM}"; then
        echo "[Rephasing-array]: FAILED samtools index for ${REGION_HAPLOTAGGED_BAM}"
        touch "${FAIL_MARKER}"
        continue
    fi

    touch "${SUCCESS_MARKER}"

    echo "[Rephasing-array]: Finished BED row ${REGION_INDEX}: ${REGION}"

done

echo "[Rephasing-array]: Task ${TASK_ID} completed all assigned regions."
EOS

    chmod +x "${ARRAY_SCRIPT}"

    echo "[Rephasing]: Submitting Slurm array with ${N_ARRAY_TASKS} total tasks"

    local ARRAY_JOB_ID
    ARRAY_JOB_ID=$(sbatch \
        --parsable \
        --job-name="rephase_${SAMPLE_CHILD}" \
        --account="${ACCOUNT}" \
        ${PARTITION:+--partition="${PARTITION}"} \
        --nodes=1 \
        --ntasks=1 \
        --array="${ARRAY_SPEC}" \
        --cpus-per-task="${REPHASE_CPUS}" \
        --mem="${REPHASE_MEM}" \
        --time="${REPHASE_TIME}" \
        --output="${LOG_DIR}/rephase_%A_%a.out" \
        --error="${LOG_DIR}/rephase_%A_%a.err" \
        "${ARRAY_SCRIPT}" "${CONFIG_FILE}" "${N_ARRAY_TASKS}")

    echo "[Rephasing]: Submitted array job: ${ARRAY_JOB_ID}"

    cat > "${SUMMARY_SCRIPT}" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="$1"
TOTAL_REGIONS="$2"

source "${CONFIG_FILE}"

SUCCESS=$(find "${LOCAL_VCF_DIR}" -name 'region_*.success' | wc -l)
FAILED=$(find "${LOCAL_VCF_DIR}" -name 'region_*.failed' | wc -l)
OBSERVED=$((SUCCESS + FAILED))
MISSING=$((TOTAL_REGIONS - OBSERVED))

echo "[Rephasing-summary]: Done"
echo "[Rephasing-summary]: Regions attempted: ${TOTAL_REGIONS}"
echo "[Rephasing-summary]: Regions successfully rephased or skipped-empty: ${SUCCESS}"
echo "[Rephasing-summary]: Regions failed: ${FAILED}"
echo "[Rephasing-summary]: Regions missing marker: ${MISSING}"
echo "[Rephasing-summary]: Local phased VCFs kept in: ${LOCAL_VCF_DIR}"

SUMMARY_FILE="${REPHASE_DIR}/rephase_summary.txt"

{
    echo "Regions attempted: ${TOTAL_REGIONS}"
    echo "Regions successfully rephased or skipped-empty: ${SUCCESS}"
    echo "Regions failed: ${FAILED}"
    echo "Regions missing marker: ${MISSING}"
    echo "Local phased VCFs: ${LOCAL_VCF_DIR}"
} > "${SUMMARY_FILE}"

echo "[Rephasing-summary]: Wrote ${SUMMARY_FILE}"

if [[ "${MISSING}" -ne 0 ]]; then
    echo "[Rephasing-summary]: WARNING: Some regions did not produce success or failure markers."
fi
EOS

    chmod +x "${SUMMARY_SCRIPT}"

    local SUMMARY_JOB_ID
    SUMMARY_JOB_ID=$(sbatch \
        --parsable \
        --dependency=afterany:${ARRAY_JOB_ID} \
        --job-name="rephase_summary_${SAMPLE_CHILD}" \
        --account="${ACCOUNT}" \
        ${PARTITION:+--partition="${PARTITION}"} \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem="2G" \
        --time="00:35:00" \
        --output="${LOG_DIR}/rephase_summary_%j.out" \
        --error="${LOG_DIR}/rephase_summary_%j.err" \
        "${SUMMARY_SCRIPT}" "${CONFIG_FILE}" "${TOTAL_REGIONS}")

    echo "[Rephasing]: Submitted summary job: ${SUMMARY_JOB_ID}"
    echo "[Rephasing]: Logs will be written to: ${LOG_DIR}"
}




# ===========================================================================
## STEP 6 : Merge rephased child VCF with unphased parent VCF, extract phased SNPs, recount mismatches
## ==========================================================================

run_rephase_merge() {
    echo "[6a-rephase] Merging parent VCF with locally rephased child VCF"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local REPHASE_DIR="${PHASED_DIR}/rephased_blocks"
    local LOCAL_VCF_DIR="${REPHASE_DIR}/local_vcfs"
    local MERGED_DIR="${REPHASE_DIR}/merged"

    mkdir -p "${MERGED_DIR}"

    if command -v module &> /dev/null; then
        module load bcftools || true
    fi

    # Part 6a depends on the locally rephased VCFs produced by Part 5
    # (rephase_blocks_with_mismatches), which writes into local_vcfs/ via a
    # Slurm array. If that directory is missing, Part 5 has not been run, was
    # skipped (no block crossed --threshold_rephase), or its array jobs have
    # not finished yet.
    if [[ ! -d "${LOCAL_VCF_DIR}" ]]; then
        echo "ERROR: ${LOCAL_VCF_DIR} does not exist."
        echo "       Part 6a needs the rephased VCFs produced by Part 5."
        echo "       Run Part 5 first:  ./main.sh --part 5 ... (same --prj-dir/--sample-*)"
        echo "       and wait for its Slurm array (rephase_${SAMPLE_CHILD}) + summary job to finish."
        echo "       Check ${REPHASE_DIR}/rephase_summary.txt for completion status."
        return 1
    fi

    ## Remove any files that's not phased or haplotagged
    find "${LOCAL_VCF_DIR}" -maxdepth 1 -type f ! -name "*phased*" ! -name "*haplotagged*" -delete

    # If Part 5 ran but produced no phased VCFs (e.g. all regions empty/failed),
    # there is nothing to merge. Fail clearly rather than silently producing
    # an empty rephase candidate set downstream.
    shopt -s nullglob
    local _have_phased=("${LOCAL_VCF_DIR}"/*phased.vcf.gz)
    shopt -u nullglob
    if [[ "${#_have_phased[@]}" -eq 0 ]]; then
        echo "ERROR: no *phased.vcf.gz files in ${LOCAL_VCF_DIR}."
        echo "       Part 5 may still be running, or no blocks exceeded --threshold_rephase,"
        echo "       or all array tasks failed. See ${REPHASE_DIR}/rephase_summary.txt"
        echo "       and ${REPHASE_DIR}/logs/ for details."
        return 1
    fi

    ## For all phased region vcfs
    shopt -s nullglob

    for phased_vcf in "${LOCAL_VCF_DIR}"/*phased.vcf.gz; do

        region_prefix=$(basename "${phased_vcf}" .phased.vcf.gz)

        merged_vcf="${MERGED_DIR}/${SAMPLE_CHILD}_${SAMPLE_PARENT}_${region_prefix}.rephased.merged.vcf.gz"

        parent_vcf="${ILLUM_DIR}/${SAMPLE_PARENT}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz"

        parent_region="${MERGED_DIR}/${region_prefix}.parent.region.vcf.gz"

        echo "Processing ${region_prefix}"

        # detect region from child phased VCF
        region=$(bcftools query -f '%CHROM\t%POS\n' "${phased_vcf}" \
            | awk '
                NR==1 {chr=$1; min=$2; max=$2}
                {
                    if($2<min) min=$2
                    if($2>max) max=$2
                }
                END {print chr ":" min "-" max}
            ')

        echo "Region wanted in the child's rephased vcf = ${region}"

        # subset parent only there
        echo "Subsetting this region from the parent's vcf" 
        bcftools view \
            -r "${region}" \
            -Oz \
            -o "${parent_region}" \
            "${parent_vcf}"

        bcftools index -f "${parent_region}"

        # merge regional files
        echo "Merging the subsetted region in the child and the parent vcfs"
        bcftools merge -m none -Oz \
            -o "${merged_vcf}" \
            "${parent_region}" \
            "${phased_vcf}"

        bcftools index -f "${merged_vcf}"

        # Remove regional parent VCF to save space
        rm -f "${parent_region}" "${parent_region}.csi" 
        echo
    done
    echo "[6a-rephase] Done merging. Merged VCFs in: ${MERGED_DIR}"
}


# ===========================================================================
## STEP 6a (Slurm): submit the remerge as ONE job (8h), then chain the rest
##   as a dependent job (afterok). Both jobs re-invoke THIS script with
##   internal stage tokens so the function logic is single-sourced.
## ==========================================================================

# Reconstruct the exact argument set needed to re-invoke main.sh inside Slurm.
# Only the flags that affect 6a/6b/6d behaviour are forwarded.
_rephase_common_args() {
    printf '%s' \
"--prj-dir ${PRJ_DIR} \
--sample-child ${SAMPLE_CHILD} \
--sample-parent ${SAMPLE_PARENT} \
--source $(basename "${SOURCE_FILE}") \
--min-rdepth ${MIN_RDEPTH} \
--max-rdepth ${MAX_RDEPTH} \
--gt-qual ${GT_QUAL} \
--nv-quantile ${NV_QUANTILE} \
--mm-diff-min ${MM_DIFF_MIN} \
--min-base-qual ${MIN_BASE_QUAL} \
--min-map-qual ${MIN_MAP_QUAL} \
--window ${WINDW} \
--alt-read-count ${ALT_READ_COUNT} \
--verbose ${VERBOSE} \
--window_rephase ${WINDOW_REPHASE} \
--threshold_rephase ${THRESHOLD_REPHASE} \
--conda-env ${CONDA_ENV_NAME}"
}

submit_rephase_6a_jobs() {
    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local REPHASE_DIR="${PHASED_DIR}/rephased_blocks"
    local LOCAL_VCF_DIR="${REPHASE_DIR}/local_vcfs"
    local LOG_DIR="${REPHASE_DIR}/logs"
    local SCRIPT_PATH
    SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

    mkdir -p "${LOG_DIR}"
    mkdir -p slurm_output

    # --- Preflight: same checks the merge loop would do, but fail fast here ---
    if [[ ! -d "${LOCAL_VCF_DIR}" ]]; then
        echo "ERROR: ${LOCAL_VCF_DIR} does not exist."
        echo "       Part 6a needs the rephased VCFs produced by Part 5."
        echo "       Run Part 5 first and let its Slurm array + summary job finish."
        return 1
    fi

    shopt -s nullglob
    local _have_phased=("${LOCAL_VCF_DIR}"/*phased.vcf.gz)
    shopt -u nullglob
    if [[ "${#_have_phased[@]}" -eq 0 ]]; then
        echo "ERROR: no *phased.vcf.gz files in ${LOCAL_VCF_DIR}."
        echo "       Part 5 may still be running, or produced no rephased regions."
        echo "       See ${REPHASE_DIR}/rephase_summary.txt"
        return 1
    fi
    echo "[6a-submit] Found ${#_have_phased[@]} rephased region VCFs to merge."

    local COMMON_ARGS
    COMMON_ARGS="$(_rephase_common_args)"

    # ---------------- Job 1: the merge (8h walltime) ----------------
    local MERGE_SLURM="${REPHASE_DIR}/run_6a_merge_${SAMPLE_CHILD}.slurm"

    cat << EOF > "${MERGE_SLURM}"
#!/bin/bash
#SBATCH -J merge6a_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH -A ${ACCOUNT}
${SBATCH_PARTITION_LINE}
#SBATCH -o ${LOG_DIR}/merge6a_%j.out
#SBATCH -e ${LOG_DIR}/merge6a_%j.err

set -euo pipefail
echo "[6a-merge-job] Host: \$(hostname)  Start: \$(date)"
bash "${SCRIPT_PATH}" ${COMMON_ARGS} --part 6a-merge
echo "[6a-merge-job] Done: \$(date)"
EOF

    chmod +x "${MERGE_SLURM}"

    local MERGE_JOBID
    MERGE_JOBID=$(sbatch --parsable "${MERGE_SLURM}")
    echo "[6a-submit] Submitted MERGE job: ${MERGE_JOBID}"

    # ---------------- Job 2: the rest (depends on merge OK) ----------------
    local REST_SLURM="${REPHASE_DIR}/run_6a_rest_${SAMPLE_CHILD}.slurm"

    cat << EOF > "${REST_SLURM}"
#!/bin/bash
#SBATCH -J rest6a_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --mem=12G
#SBATCH -A ${ACCOUNT}
${SBATCH_PARTITION_LINE}
#SBATCH -o ${LOG_DIR}/rest6a_%j.out
#SBATCH -e ${LOG_DIR}/rest6a_%j.err

set -euo pipefail
echo "[6a-rest-job] Host: \$(hostname)  Start: \$(date)"
bash "${SCRIPT_PATH}" ${COMMON_ARGS} --part 6a-rest
echo "[6a-rest-job] Done: \$(date)"
EOF

    chmod +x "${REST_SLURM}"

    local REST_JOBID
    REST_JOBID=$(sbatch --parsable --dependency=afterok:${MERGE_JOBID} "${REST_SLURM}")
    echo "[6a-submit] Submitted REST job: ${REST_JOBID} (runs after ${MERGE_JOBID} succeeds)"

    echo "[6a-submit] Merge script : ${MERGE_SLURM}"
    echo "[6a-submit] Rest script  : ${REST_SLURM}"
    echo "[6a-submit] Logs         : ${LOG_DIR}"
    echo "[6a-submit] Monitor with : squeue -u \$USER"
}

extract_rephase_phased_snp() {
    echo "[6a-rephase] Extracting phased SNPs from ALL rephased merged VCFs"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local REPHASE_DIR="${PHASED_DIR}/rephased_blocks"
    local MERGED_DIR="${REPHASE_DIR}/merged"
    local HT_DIR="${PHASED_DIR}/HTblocks"

    local OUT_TSV="${HT_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps${REPHASE_SUFFIX}.tsv"
    local FIXED_TSV="${HT_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps${REPHASE_SUFFIX}_fixed.tsv"

    mkdir -p "${HT_DIR}"

    # clear previous file
    : > "${OUT_TSV}"

    shopt -s nullglob
    local found=0

    for merged_vcf in "${MERGED_DIR}"/*.vcf.gz; do
        found=1
        echo "[6a-rephase] Processing $(basename "${merged_vcf}")"

        bcftools view "${merged_vcf}" \
            -i 'TYPE="snp" && N_ALT=1 && GT[0]!="mis" && GT[1]!="mis" && FMT/PS[1]!="."' \
            -Ou \
        | bcftools query \
            -s "${SAMPLE_PARENT},${SAMPLE_CHILD}" \
            -f '%CHROM,%POS,[%PS,][%GT,][%DP,][%GQ,][%AD{0},][%AD{1},]\n' \
            >> "${OUT_TSV}"
    done
    shopt -u nullglob

    if [[ "${found}" -eq 0 ]]; then
        echo "[6a-rephase] No merged VCF files found in ${MERGED_DIR}"
        return 1
    fi

    if [[ ! -s "${OUT_TSV}" ]]; then
        echo "[6a-rephase] No phased SNPs found"
        return 1
    fi

    # Sort + unique by chromosome/position
    sort -t',' -k1,1 -k2,2n -u "${OUT_TSV}" > "${OUT_TSV}.tmp"
    mv "${OUT_TSV}.tmp" "${OUT_TSV}"

    conda activate "${CONDA_ENV_NAME}"

    python "${WORKING_DIR}/src/fix_PhaseSet.py" \
        "${OUT_TSV}" \
        "${FIXED_TSV}"

    mv "${FIXED_TSV}" "${OUT_TSV}"

    echo "[6a-rephase] Final concatenated output:"
    echo "${OUT_TSV}"
}

count_rephase_mismatches() {
    echo "[6a-rephase] Counting mismatches using rephased blocks"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local HT_DIR="${PHASED_DIR}/HTblocks"
    local MISMATCH_ANALYSIS_DIR="${PHASED_DIR}/mismatch_analysis"

    
    conda activate "${CONDA_ENV_NAME}"

    python "${WORKING_DIR}/src/count_rephase_mismatches.py" \
    "${HT_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps${REPHASE_SUFFIX}.tsv" \
    "${MISMATCH_ANALYSIS_DIR}" \
    "${MIN_RDEPTH}" \
    "${MAX_RDEPTH}" \
    "${SAMPLE_PARENT}" \
    "${SAMPLE_CHILD}" \
    "${GT_QUAL}" \
    "${NV_QUANTILE}" \
    "${MM_DIFF_MIN}" 


    echo "[6a-rephase] See output in file ${MISMATCH_ANALYSIS_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_mismatch${REPHASE_SUFFIX}.tsv"
    echo "[6a-rephase] Done counting mismatches with rephased blocks"
}


filter_rephase_dnm_candidates() {
    echo "[6a-rephase] Filtering rephased DNM candidates"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local MISMATCH_ANALYSIS_DIR="${PHASED_DIR}/mismatch_analysis"
    local REPHASED_MERGED_DIR="${PHASED_DIR}/rephased_blocks/merged"

    local IN_BED="${MISMATCH_ANALYSIS_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc${REPHASE_SUFFIX}.bed"
    local OUT_TSV="${MISMATCH_ANALYSIS_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc${REPHASE_SUFFIX}.tsv"
    local OUT_BED="${MISMATCH_ANALYSIS_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc${REPHASE_SUFFIX}.filt.bed"

    # Clear outputs before appending
    : > "${OUT_TSV}"
    : > "${OUT_BED}"

    # Hardened filter expression
    # [0] = parent (must be hom-ref, zero alt reads, passing GQ+DP)
    # [1] = child  (must be phased het, passing GQ+DP, alt reads present)
    local filter_expr="
    (FMT/GT[1]==\"0|1\" || FMT/GT[1]==\"1|0\") &&
    FMT/GT[0]==\"0/0\" &&
    FMT/GT[0]!=\"1/1\" &&
    FMT/GT[0]!=\"0/1\" &&
    FMT/GT[0]!=\"1/0\" &&
    FMT/DP[0]>=${MIN_RDEPTH} && FMT/DP[0]<=${MAX_RDEPTH} &&
    FMT/DP[1]>=${MIN_RDEPTH} && FMT/DP[1]<=${MAX_RDEPTH} &&
    FMT/GQ[0]>=${GT_QUAL} &&
    FMT/GQ[1]>=${GT_QUAL} &&
    FMT/AD[1:1]>5 &&
    FMT/AD[1:0]>5 &&
    FMT/AD[0:1]==0
    "

    shopt -s nullglob
    local found=0
    local header_written=0

    for merged_vcf in "${REPHASED_MERGED_DIR}"/*.vcf.gz; do
        found=1
        echo "[6a-rephase] Processing $(basename "${merged_vcf}")"

        local tmp_vcf="${OUT_TSV%.tsv}.temp.vcf.gz"

        bcftools view \
            -R "${IN_BED}" \
            -m2 -M2 -v snps \
            -Oz \
            -o "${tmp_vcf}" \
            "${merged_vcf}"

        bcftools index -f "${tmp_vcf}"

        # Verify sample order: parent must be column 0
        local col0
        col0=$(bcftools query -l "${tmp_vcf}" | sed -n '1p')
        if [[ "${col0}" != "${SAMPLE_PARENT}" ]]; then
            echo "ERROR: expected ${SAMPLE_PARENT} as sample[0] in ${merged_vcf}, got ${col0}. Skipping."
            rm -f "${tmp_vcf}" "${tmp_vcf}.csi"
            continue
        fi

        # Write the bcftools -H header exactly ONCE (from the first valid VCF),
        # then append data rows WITHOUT -H so the header is not repeated per region.
        if [[ "${header_written}" -eq 0 ]]; then
            bcftools view "${tmp_vcf}" -i "${filter_expr}" \
            | bcftools query -H \
                -f '%CHROM\t%POS\t[%GT\t][%DP\t][%GQ\t][%AD\t]\n' \
            | head -n 1 \
            > "${OUT_TSV}"
            header_written=1
        fi

        # Data rows only (no -H): bcftools query without -H emits no header line.
        bcftools view "${tmp_vcf}" -i "${filter_expr}" \
        | bcftools query \
            -f '%CHROM\t%POS\t[%GT\t][%DP\t][%GQ\t][%AD\t]\n' \
            >> "${OUT_TSV}"

        bcftools view "${tmp_vcf}" -i "${filter_expr}" \
        | bcftools query -f '%CHROM\t%POS\n' \
        | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' \
            >> "${OUT_BED}"

        rm -f "${tmp_vcf}" "${tmp_vcf}.csi"
    done
    shopt -u nullglob

    if [[ "${found}" -eq 0 ]]; then
        echo "[6a-rephase] WARNING: No merged VCFs found in ${REPHASED_MERGED_DIR}"
        return 1
    fi

    # De-duplicate and sort the BED
    if [[ -s "${OUT_BED}" ]]; then
        sort -k1,1 -k2,2n -u "${OUT_BED}" > "${OUT_BED}.tmp"
        mv "${OUT_BED}.tmp" "${OUT_BED}"
    fi

    echo "[6a-rephase] Input BED  : ${IN_BED}"
    echo "[6a-rephase] Output TSV : ${OUT_TSV}"
    echo "[6a-rephase] Output BED : ${OUT_BED}"
    echo "[6a-rephase] Candidate count: $(wc -l < "${OUT_BED}")"
}


validate_rephase_dnmc_with_parent_bam() {
    echo "[6a-parent-check] Validating rephased DNM candidates against parent Illumina BAM"

    conda activate "${CONDA_ENV_NAME}"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local MISMATCH_ANALYSIS_DIR="${PHASED_DIR}/mismatch_analysis"
    local MERGED_PHASED_VCF="${PHASED_DIR}/merged"

    # The rephase candidate BED that filter_rephase_dnm_candidates produced and
    # that the HiFi validation will consume next. We filter it IN PLACE, exactly
    # as the original 2b validate_dnmc_with_parent_bam does for set A.
    local IN_BED="${MISMATCH_ANALYSIS_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc${REPHASE_SUFFIX}.filt.bed"
    local OUT_BED="${MISMATCH_ANALYSIS_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc${REPHASE_SUFFIX}.parentBAM.bed"

    # parent_readcheck.py looks up REF/ALT per site via vcf.fetch(), so the VCF
    # must contain every candidate position. The per-region rephase merged VCFs
    # do NOT collectively guarantee that in one file, but the whole-genome 2b
    # merged VCF does, and the parent's REF/ALT at any site are identical there
    # (rephasing only changed the CHILD's phasing, not the parent record). So we
    # use the canonical 2b merged VCF.
    local MERGED_VCF="${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz"
    local PARENT_BAM="${ORIG_BAM_PARENT_PATH}"
    local MAX_PARENT_ALT=1

    if [[ ! -s "${IN_BED}" ]]; then
        echo "[6a-parent-check] No candidates in ${IN_BED}; nothing to do"
        return 0
    fi

    if [[ ! -f "${MERGED_VCF}" ]]; then
        echo "[6a-parent-check] WARNING: 2b merged VCF not found at ${MERGED_VCF}"
        echo "[6a-parent-check] Skipping parent-BAM validation (run Parts 2b-4 first)."
        return 0
    fi

    if [[ ! -f "${PARENT_BAM}" ]]; then
        echo "[6a-parent-check] WARNING: parent BAM not found at ${PARENT_BAM}"
        echo "[6a-parent-check] Skipping parent-BAM validation (results may include FN-parent DNMs)"
        return 0
    fi

    # Make sure the parent BAM is indexed
    if [[ ! -f "${PARENT_BAM}.bai" && ! -f "${PARENT_BAM%.bam}.bai" ]]; then
        echo "[6a-parent-check] Indexing parent BAM"
        samtools index "${PARENT_BAM}"
    fi

    local N_BEFORE
    N_BEFORE=$(grep -cv '^#' "${IN_BED}" || true)

    python "${WORKING_DIR}/src/parent_readcheck.py" \
        "${SAMPLE_PARENT}" \
        "${PARENT_BAM}" \
        "${MERGED_VCF}" \
        "${IN_BED}" \
        "${OUT_BED}" \
        "${MIN_BASE_QUAL}" \
        "${MIN_MAP_QUAL}" \
        "${MAX_PARENT_ALT}" \
        "${MISMATCH_ANALYSIS_DIR}"

    # Swap filtered BED into the canonical location (so HiFi validation sees it)
    mv "${OUT_BED}" "${IN_BED}"

    local N_AFTER
    N_AFTER=$(grep -cv '^#' "${IN_BED}" || true)

    echo "[6a-parent-check] Candidates before: ${N_BEFORE}"
    echo "[6a-parent-check] Candidates after : ${N_AFTER}"
    echo "[6a-parent-check] Removed by parent-BAM evidence: $((N_BEFORE - N_AFTER))"
}


validate_rephase_dnmc_with_hifi_reads() {
    echo "[6a-rephase] Validating rephased DNMC candidates with HiFi reads"

    conda activate "${CONDA_ENV_NAME}"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local REPHASE_ANALYSIS_DIR="${PHASED_DIR}/mismatch_analysis"

    # Use the bcftools-filtered output rather than the raw Python output
    local DNM_BED="${REPHASE_ANALYSIS_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc${REPHASE_SUFFIX}.filt.bed"

    # Distinct label avoids output-file collisions with the set-A 'dnmc' run
    local LR_LABEL="dnmc${REPHASE_SUFFIX}"

    # dnmc_readcheck.py now writes directly into REPHASE_ANALYSIS_DIR (arg 11),
    # so there is no CWD-to-dir mv afterward.
    python "${WORKING_DIR}/src/dnmc_readcheck.py" \
        "${SAMPLE_CHILD}" \
        "${HP_BAM_PATH}" \
        "${DNM_BED}" \
        "${MIN_BASE_QUAL}" \
        "${MIN_MAP_QUAL}" \
        "${WINDW}" \
        "${ALT_READ_COUNT}" \
        "T" \
        "${LR_LABEL}" \
        "${TOTAL_READ_COUNT_MIN}" \
        "${REPHASE_ANALYSIS_DIR}"

    # ----------------------------------------------------------------------
    # Make the LR-validated set the CANONICAL set-B candidate list, mirroring
    # how set A treats LR validation as a hard filter on candidates.
    #
    # dnmc_readcheck.py wrote ${CHILD}_LR_validated_${LR_LABEL}.bed directly
    # into REPHASE_ANALYSIS_DIR. 
    ## Here, we intersect the bcftools-filtered candidate TSV against the LR-validated positions so that
    #   ${PARENT}_${CHILD}_dnmc_rephase.tsv
    # contains ONLY LR-passing candidates (header preserved).
    # ----------------------------------------------------------------------
    local LR_VALID_BED="${REPHASE_ANALYSIS_DIR}/${SAMPLE_CHILD}_LR_validated_${LR_LABEL}.bed"
    local FILT_TSV="${REPHASE_ANALYSIS_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc${REPHASE_SUFFIX}.tsv"

    if [[ -s "${LR_VALID_BED}" && -s "${FILT_TSV}" ]]; then
        echo "[6a-rephase] Intersecting bcftools-filtered candidates with LR-validated positions"

        local TMP_TSV="${FILT_TSV%.tsv}.lrvalidated.tsv"

        awk '
            FNR==NR { v[$1":"$3]=1; next }
            FNR==1 { print; next }
            $1 ~ /^#/ { next }
            {
                key=$1":"$2
                if (key in v) print
            }
        ' "${LR_VALID_BED}" "${FILT_TSV}" > "${TMP_TSV}"

        mv "${TMP_TSV}" "${FILT_TSV}"

        echo "[6a-rephase] Canonical LR-validated set-B candidates: ${FILT_TSV}"
        echo "[6a-rephase] Set-B candidate count (LR-validated): $(awk 'NR>1 && $1 !~ /^#/ && NF>0' "${FILT_TSV}" | wc -l)"
    else
        echo "[6a-rephase] No LR-validated candidates: setting set B to ZERO."
        echo "[6a-rephase]   LR BED : ${LR_VALID_BED}"
        echo "[6a-rephase]   TSV    : ${FILT_TSV}"

        if [[ -s "${FILT_TSV}" ]]; then
            head -n 1 "${FILT_TSV}" > "${FILT_TSV}.tmp"
            mv "${FILT_TSV}.tmp" "${FILT_TSV}"
        else
            printf '# [1]CHROM\t[2]POS\n' > "${FILT_TSV}"
        fi
        echo "[6a-rephase] Set-B candidate count (LR-validated): 0"
    fi
}

# ===========================================================================
## STEP 6b : INDEPENDENT callable genome over the rephased blocks (gives C_B)
##           Mirrors Part 4 (calculate_callable_genome + REfilter + REvalidate)
##           but points at the rephased PS file / rephased merged VCFs, writing
##           to a SEPARATE denum dir so set A's callable_genome.txt is untouched.
## ==========================================================================

calculate_callable_genome_rephase() {
    echo "[6b-rephase] Computing INDEPENDENT callable genome over rephased blocks"

    conda activate "${CONDA_ENV_NAME}"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local HT_DIR="${PHASED_DIR}/HTblocks"
    local DENUM_DIR="${PHASED_DIR}/mismatch_analysis${REPHASE_SUFFIX}/denum_calcul"

    # Rephased PS file produced by extract_rephase_phased_snp (Part 6a)
    local REPHASE_PS="${HT_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps${REPHASE_SUFFIX}.tsv"

    if [[ ! -s "${REPHASE_PS}" ]]; then
        echo "ERROR: rephased PS file not found or empty: ${REPHASE_PS}"
        return 1
    fi

    mkdir -p "${DENUM_DIR}"

    EXCLUDE_CHROMS="${EXCLUDE_CHROMS}" \
    # Sample callable SNPs from the rephased blocks (same args as Part 4)
    python "${WORKING_DIR}/src/callable_genome.py" \
        "${REPHASE_PS}" \
        "${DENUM_DIR}" \
        "${MIN_RDEPTH}" \
        "${MAX_RDEPTH}" \
        "${SAMPLE_PARENT}" \
        "${SAMPLE_CHILD}" \
        "${GT_QUAL}" \
        "${NV_QUANTILE}" \
        "${MM_DIFF_MIN}" \
        0 "${NOTRECOUNT}"

    echo "[6b-rephase] Independent callable-genome outputs in: ${DENUM_DIR}"
}


REfilter_dnm_candidates_rephase() {
    echo "[6b-rephase] Re-filtering sampled SNPs over rephased merged VCFs (for C_B)"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local DENUM_DIR="${PHASED_DIR}/mismatch_analysis${REPHASE_SUFFIX}/denum_calcul"
    local REPHASED_MERGED_DIR="${PHASED_DIR}/rephased_blocks/merged"

    local BED="${DENUM_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.bed"
    local OUT_TSV="${DENUM_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.tsv"
    local OUT_BED="${DENUM_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.filt.bed"

    if [[ ! -s "${BED}" ]]; then
        echo "ERROR: het-candidate BED not found or empty: ${BED}"
        return 1
    fi

    : > "${OUT_TSV}"
    : > "${OUT_BED}"

    local filter_expr="
    (FMT/GT[1]==\"0|1\" || FMT/GT[1]==\"1|0\") &&
    (FMT/GT[0]==\"0/1\" || FMT/GT[0]==\"1/0\") &&
    FMT/DP>=${MIN_RDEPTH} && FMT/DP<=${MAX_RDEPTH} &&
    FMT/GQ>=${GT_QUAL} &&
    FMT/AD[1:1]>5 &&
    FMT/AD[1:0]>5
    "

    shopt -s nullglob
    local found=0
    local header_written=0
    for merged_vcf in "${REPHASED_MERGED_DIR}"/*.vcf.gz; do
        found=1
        local tmp_vcf="${OUT_TSV%.tsv}.temp.vcf.gz"

        bcftools view -R "${BED}" -m2 -M2 -v snps -Oz -o "${tmp_vcf}" "${merged_vcf}"
        bcftools index -f "${tmp_vcf}"

        local col0
        col0=$(bcftools query -l "${tmp_vcf}" | sed -n '1p')
        if [[ "${col0}" != "${SAMPLE_PARENT}" ]]; then
            echo "WARN: ${SAMPLE_PARENT} not sample[0] in $(basename "${merged_vcf}"); skipping."
            rm -f "${tmp_vcf}" "${tmp_vcf}.csi"; continue
        fi

        # Header exactly once (from first valid VCF); data rows without -H after.
        if [[ "${header_written}" -eq 0 ]]; then
            bcftools view "${tmp_vcf}" -i "${filter_expr}" \
            | bcftools query -H -f '%CHROM\t%POS\t[%GT\t][%DP\t][%GQ\t][%AD\t]\n' \
            | head -n 1 \
            > "${OUT_TSV}"
            header_written=1
        fi

        bcftools view "${tmp_vcf}" -i "${filter_expr}" \
        | bcftools query -f '%CHROM\t%POS\t[%GT\t][%DP\t][%GQ\t][%AD\t]\n' \
            >> "${OUT_TSV}"

        bcftools view "${tmp_vcf}" -i "${filter_expr}" \
        | bcftools query -f '%CHROM\t%POS\n' \
        | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' \
            >> "${OUT_BED}"

        rm -f "${tmp_vcf}" "${tmp_vcf}.csi"
    done
    shopt -u nullglob

    if [[ "${found}" -eq 0 ]]; then
        echo "ERROR: no merged VCFs in ${REPHASED_MERGED_DIR}"
        return 1
    fi

    if [[ -s "${OUT_BED}" ]]; then
        sort -k1,1 -k2,2n -u "${OUT_BED}" > "${OUT_BED}.tmp"
        mv "${OUT_BED}.tmp" "${OUT_BED}"
    fi

    # canonicalize name expected by LR validation (overwrite the sampled BED)
    mv "${OUT_BED}" "${BED}"
    echo "[6b-rephase] Qualified-SNP BED (pre-LR): ${BED}"
}



REvalidate_dnmc_with_long_reads_rephase() {
    echo "[6b-rephase] LR-validating sampled SNPs over rephased blocks (for C_B)"
    conda activate "${CONDA_ENV_NAME}"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local DENUM_DIR="${PHASED_DIR}/mismatch_analysis${REPHASE_SUFFIX}/denum_calcul"
    local DNM_BED="${DENUM_DIR}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.bed"

    # Distinct label avoids collisions with Part 4's 'hetc' outputs
    local LR_LABEL="hetc${REPHASE_SUFFIX}"

    # dnmc_readcheck.py now writes directly into DENUM_DIR (arg 11); no mv after.
    python "${WORKING_DIR}/src/dnmc_readcheck.py" \
        "${SAMPLE_CHILD}" \
        "${HP_BAM_PATH}" \
        "${DNM_BED}" \
        "${MIN_BASE_QUAL}" \
        "${MIN_MAP_QUAL}" \
        "${WINDW}" \
        "${ALT_READ_COUNT}" \
        "F" \
        "${LR_LABEL}" \
        "${TOTAL_READ_COUNT_MIN}" \
        "${DENUM_DIR}"

    echo "[6b-rephase] Independent LR-validated qualified SNPs:"
    echo "             ${DENUM_DIR}/${SAMPLE_CHILD}_LR_validated_${LR_LABEL}.bed"
}

# ===========================================================================
## STEP 6d (zero path): No LR-validated rephase DNMs.
##   Set B contributes 0 candidates, so the candidate-weighted rate collapses
##   to the original set-A rate. We do NOT compute C_B (pointless when d_B=0).
##   The merged file == set A; the NEW-only file is empty (header only).
## ==========================================================================

final_summary_rephase_zero() {
    echo "[6d-rephase] Set B = 0 validated DNMs; final rate = original set-A rate"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"

    # --- Set A (previous) ---
    local A_DNMC="${PHASED_DIR}/final_dnmc_${SAMPLE_CHILD}-from-${SAMPLE_PARENT}.tsv"
    local A_DENUM="${PHASED_DIR}/mismatch_analysis/denum_calcul"
    local A_CALLABLE="${A_DENUM}/callable_genome.txt"
    local A_QUALBED="${A_DENUM}/${SAMPLE_CHILD}_LR_validated_hetc.bed"

    local MERGED_DNMC="${PHASED_DIR}/final_dnmc_${SAMPLE_CHILD}-from-${SAMPLE_PARENT}_with_rephase.tsv"
    local B_NEW_DNMC="${PHASED_DIR}/final_dnmc_${SAMPLE_CHILD}-from-${SAMPLE_PARENT}_rephase_NEWonly.tsv"
    local SUMMARY_FILE="${PHASED_DIR}/final_summary_with_rephase.txt"

    if [[ ! -f "${A_DNMC}" ]]; then
        echo "ERROR: set-A final DNM file missing: ${A_DNMC}"
        return 1
    fi
    if [[ ! -f "${A_CALLABLE}" ]]; then
        echo "ERROR: set-A callable_genome.txt missing: ${A_CALLABLE}"
        return 1
    fi

    # Merged == set A (nothing new to add); NEW-only is empty (header only).
    cp -f "${A_DNMC}" "${MERGED_DNMC}"
    head -n 1 "${A_DNMC}" > "${B_NEW_DNMC}"

    # ---- Set A counts and rate (identical to Part 4's final_summary) ----
    local dA qual_A sampled_A accessible_A
    dA=$(awk 'NR>1 && $1 !~ /^#/ && NF>0' "${A_DNMC}" | wc -l)
    qual_A=$(wc -l < "${A_QUALBED}")
    sampled_A=$(grep "total_sampled_snps"    "${A_CALLABLE}" | awk '{print $2}')
    accessible_A=$(grep "total_callable_bases" "${A_CALLABLE}" | awk '{print $2}')

    local C_A r_A
    read -r C_A r_A <<<"$(awk \
        -v dA="${dA}" -v qA="${qual_A}" -v sA="${sampled_A:-0}" -v cA="${accessible_A:-0}" '
        BEGIN{
            CA = (sA>0)? (qA/sA)*cA : 0
            rA = (CA>0)? dA/CA : 0
            printf "%.6e %.6e", CA, rA
        }')"

    {
        echo "================ Final results: previous + rephase (set B = 0) ============"
        echo "Independent file (set A) : ${A_DNMC}"
        echo "Set B                    : 0 LR-validated DNMs (no contribution)"
        echo "Set B new-only file      : ${B_NEW_DNMC} (empty)"
        echo "Merged file (A union B)  : ${MERGED_DNMC} (== set A)"
        echo "-------------------------------------------------------------------"
        echo "Set A : d_A=${dA}  sampled=${sampled_A}  qualified=${qual_A}  accessible=${accessible_A}"
        printf "Set A : C_A=%.6e  r_A=%.6e\n" "${C_A}" "${r_A}"
        echo "Set B : d_B=0  (skipped callable-genome-rephase)"
        echo "-------------------------------------------------------------------"
        echo "Merged unique DNM count  : ${dA}"
        printf "Final mutation rate (= set-A rate) = %.6e\n" "${r_A}"
        echo "==========================================================================="
    } | tee "${SUMMARY_FILE}"

    echo
    echo "[6d-rephase] Summary written to: ${SUMMARY_FILE}"
}


# ===========================================================================
## STEP 6d : Independent files + merged file + candidate-weighted mutation rate
##   Set A = previous results (steps 3b/4), untouched.
##   Set B = rephase results (this part), independent callable genome.
##   d_B    = ONLY part-6 candidates NOT already present in set A.
##   r_w    = (d_A*r_A + d_B*r_B) / (d_A + d_B)   [candidate-weighted]
## ==========================================================================

final_summary_with_rephase() {
    echo "[6d-rephase] Building merged DNM set and candidate-weighted mutation rate"

    local PHASED_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"

    # --- Set A (previous) ---
    local A_DNMC="${PHASED_DIR}/final_dnmc_${SAMPLE_CHILD}-from-${SAMPLE_PARENT}.tsv"
    local A_DENUM="${PHASED_DIR}/mismatch_analysis/denum_calcul"
    local A_CALLABLE="${A_DENUM}/callable_genome.txt"
    local A_QUALBED="${A_DENUM}/${SAMPLE_CHILD}_LR_validated_hetc.bed"

    # --- Set B (rephase) : LR-validated canonical candidates from Part 6a ---
    local B_DNMC="${PHASED_DIR}/mismatch_analysis/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc${REPHASE_SUFFIX}.tsv"
    local B_DENUM="${PHASED_DIR}/mismatch_analysis${REPHASE_SUFFIX}/denum_calcul"
    local B_CALLABLE="${B_DENUM}/callable_genome.txt"
    local B_QUALBED="${B_DENUM}/${SAMPLE_CHILD}_LR_validated_hetc${REPHASE_SUFFIX}.bed"

    # --- Independent + merged candidate files ---
    local MERGED_DNMC="${PHASED_DIR}/final_dnmc_${SAMPLE_CHILD}-from-${SAMPLE_PARENT}_with_rephase.tsv"
    local B_NEW_DNMC="${PHASED_DIR}/final_dnmc_${SAMPLE_CHILD}-from-${SAMPLE_PARENT}_rephase_NEWonly.tsv"
    local SUMMARY_FILE="${PHASED_DIR}/final_summary_with_rephase.txt"

    local f
    for f in "${A_DNMC}" "${B_DNMC}" "${A_CALLABLE}" "${B_CALLABLE}"; do
        if [[ ! -f "${f}" ]]; then
            echo "ERROR: required file missing: ${f}"
            return 1
        fi
    done

    # ---- Merged (A union B), dedup by chrom:pos, header from A ----
    # Skip any stray '#'-prefixed header lines (defensive) and blank lines.
    {
        head -n 1 "${A_DNMC}"
        tail -n +2 "${A_DNMC}"
        tail -n +2 "${B_DNMC}"
    } | awk '
        NR==1 { print; next }
        $1 ~ /^#/ { next }
        NF==0 { next }
        { k=$1":"$2; if(!(k in s)){ s[k]=1; print } }
    ' > "${MERGED_DNMC}"

    # ---- Set B "new only" = B candidates whose chrom:pos absent from A ----
    {
        head -n 1 "${B_DNMC}"
        awk '
            NR==FNR { if (FNR>1 && $1 !~ /^#/ && NF>0) a[$1":"$2]=1; next }
            FNR==1  { next }
            $1 ~ /^#/ { next }
            NF==0 { next }
            { if (!(($1":"$2) in a)) print }
        ' "${A_DNMC}" "${B_DNMC}"
    } > "${B_NEW_DNMC}"

    # ---- Counts (count real data rows, not wc-l-minus-1) ----
    local dA dB_new
    dA=$(awk 'NR>1 && $1 !~ /^#/ && NF>0' "${A_DNMC}" | wc -l)
    dB_new=$(awk 'NR>1 && $1 !~ /^#/ && NF>0' "${B_NEW_DNMC}" | wc -l)
    [[ "${dA}"     -lt 0 ]] && dA=0
    [[ "${dB_new}" -lt 0 ]] && dB_new=0

    # ---- Callable genome A (C_A) ----
    local sampled_A accessible_A qual_A
    sampled_A=$(grep "total_sampled_snps"    "${A_CALLABLE}" | awk '{print $2}')
    accessible_A=$(grep "total_callable_bases" "${A_CALLABLE}" | awk '{print $2}')
    qual_A=$(wc -l < "${A_QUALBED}")

    # ---- Callable genome B (C_B), independent ----
    local sampled_B accessible_B qual_B
    sampled_B=$(grep "total_sampled_snps"    "${B_CALLABLE}" | awk '{print $2}')
    accessible_B=$(grep "total_callable_bases" "${B_CALLABLE}" | awk '{print $2}')
    qual_B=$(wc -l < "${B_QUALBED}")

    # ---- Rates and candidate-weighted average ----
    # C_A = (qual_A/sampled_A)*accessible_A ;  r_A = dA / C_A
    # C_B = (qual_B/sampled_B)*accessible_B ;  r_B = dB_new / C_B
    # r_w = (dA*r_A + dB_new*r_B) / (dA + dB_new)
    local C_A r_A C_B r_B r_w
    read -r C_A r_A C_B r_B r_w <<<"$(awk \
        -v dA="${dA}" -v qA="${qual_A}" -v sA="${sampled_A:-0}" -v cA="${accessible_A:-0}" \
        -v dB="${dB_new}" -v qB="${qual_B}" -v sB="${sampled_B:-0}" -v cB="${accessible_B:-0}" '
        BEGIN{
            CA = (sA>0)? (qA/sA)*cA : 0
            CB = (sB>0)? (qB/sB)*cB : 0
            rA = (CA>0)? dA/CA : 0
            rB = (CB>0)? dB/CB : 0
            w  = ((dA+dB)>0)? (dA*rA + dB*rB)/(dA+dB) : 0
            printf "%.6e %.6e %.6e %.6e %.6e", CA, rA, CB, rB, w
        }')"

    {
        echo "================ Final results: previous + rephase ================"
        echo "Independent file (set A) : ${A_DNMC}"
        echo "Independent file (set B) : ${B_DNMC}"
        echo "Set B new-only file      : ${B_NEW_DNMC}"
        echo "Merged file (A union B)  : ${MERGED_DNMC}"
        echo "-------------------------------------------------------------------"
        echo "Set A : d_A=${dA}  sampled=${sampled_A}  qualified=${qual_A}  accessible=${accessible_A}"
        printf "Set A : C_A=%.6e  r_A=%.6e\n" "${C_A}" "${r_A}"
        echo "Set B : d_B(new)=${dB_new}  sampled=${sampled_B}  qualified=${qual_B}  accessible=${accessible_B}"
        printf "Set B : C_B=%.6e  r_B=%.6e\n" "${C_B}" "${r_B}"
        echo "-------------------------------------------------------------------"
        echo "Merged unique DNM count  : $(awk 'NR>1 && $1 !~ /^#/ && NF>0' "${MERGED_DNMC}" | wc -l)"
        printf "Candidate-weighted mutation rate = (d_A*r_A + d_B*r_B)/(d_A+d_B) = %.6e\n" "${r_w}"
        echo "==================================================================="
    } | tee "${SUMMARY_FILE}"

    echo
    echo "[6d-rephase] Summary written to: ${SUMMARY_FILE}"
}




##################################################
# THE MAIN POINT!!!!
##################################################


main() {

    load_modules  # always

    for PART in "${PARTS[@]}"; do
        ###########################################
        # PART 0 - Print out summary of the config
        ###########################################
        if [[ "$PART" == "0" ]]; then
            echo "========== PART 0: Checking configuration =========="
            echo "Make sure that you have the following dependencies"
            cat $DEPENDENCIES_FILE
        fi

        ############################################
        # PART 1A — Data download
        ############################################
        if [[ "$PART" == "1a" ]]; then
            echo "========== PART 1A: Downloading data =========="
            DOWNLOAD_JOBID=$(download_data_job)
            echo "Download job submitted: ${DOWNLOAD_JOBID}"
        fi

        ############################################
        # PART 1B — BAM preprocessing
        ############################################
        if [[ "$PART" == "1b" ]]; then
            echo "========== PART 1B: BAM preprocessing =========="
            PREPROCESS_BAM_JOBID=$(generate_bam_preprocessing_job)
            echo "BAM preprocessing (check HP tag) job submitted: ${PREPROCESS_BAM_JOBID}"
        fi

        ############################################
        # PART 1C — VCF preprocessing
        ############################################
        if [[ "$PART" == "1c" ]]; then
            echo "========== PART 1C: VCF preprocessing =========="
            PREPROCESS_VCF_JOBID=$(generate_vcf_preprocessing_job)
            echo "VCF preprocessing job submitted: ${PREPROCESS_VCF_JOBID}"
        fi


        ############################################
        # PART 2 — Initial Phasing
        ############################################
        if [[ "$PART" == "2" ]]; then
            echo "========== PART 2: Initial phasing =========="
            generate_phasing_job  # Phasing SLURM JOB, Usually take 4 hours
        
        fi


        ############################################
        # PART 2B — First DNM detection
        ############################################
        if [[ "$PART" == "2b" ]]; then
            echo "========== PART 2B: Merge + first DNM detection =========="

            # Gate: Part 2's output must actually be phased (FORMAT/PS present)
            # before we touch it. Works regardless of which phasing program
            # produced it (whatshap, --external-phased-vcf, or hand-placed - hopefully
            echo "Checking that Part 2 produced a properly phased VCF (FORMAT/PS present)"
            if ! check_phased_vcf_has_ps; then
                echo "[2b] Aborting: phased VCF check failed."
                exit 1
            fi

            # merge phased vcf of the child to unphased vcf of parent
            echo "Re-merging parent and child VCFs for DNM detection"
            merge_unphased-parent_phased-child_vcfs

            # Extract only positions where there's a Phase Set (PS) associated
            echo "Extracting phased SNPs with PS tags for DNM candidate detection"
            extract_phased_snp

            # Accounting on PS blocks
            echo "Counting shared alleles per PS block to identify mismatch blocks"
            count_shared_alleles_per_PS_block

            # Create a list of dnm candidates
            echo "Applying filters to identify de novo mutation candidates"
            filter_dnm_candidates

            # Reject candidates where parent's Illumina BAM has ALT reads
            echo "Validating DNM candidates against parent Illumina BAM"
            validate_dnmc_with_parent_bam
     

            # Additional filters with LR suppport
            echo "Validating DNM candidates with long read support"
            validate_dnmc_with_long_reads

        fi



        # ---- inside main()'s for-loop ----

        # PART 3 branch: self-skip if 0 candidates
        if [[ "$PART" == "3" ]]; then
            echo "========== PART 3: Local rephasing =========="
            if dnm_candidates_are_empty; then
                echo "[3] 0 LR-validated DNM candidates from Part 2b; nothing to rephase. Skipping."
            else
                regenerate_phasing_job
            fi
        fi

        # PART 3-inline branch (used by the chain): self-skip if 0 candidates
        if [[ "$PART" == "3-inline" ]]; then
            echo "========== PART 3-INLINE: local phasing (in-process, for chain) =========="
            if dnm_candidates_are_empty; then
                echo "[3-inline] 0 LR-validated DNM candidates; nothing to rephase. Skipping cleanly."
            else
                _run_local_phasing_body
            fi
        fi

        # PART 3B branch: self-skip if 0 candidates
        if [[ "$PART" == "3b" ]]; then
            echo "========== PART 3B: Refined DNM detection =========="
            if dnm_candidates_are_empty; then
                echo "[3b] 0 LR-validated DNM candidates; no local phasing was done. Skipping refinement."
            else
                remerge_unphased-parent_phased-child_vcfs
                reextract_phased_snp
                recount_shared_alleles_per_PS_block
                clean_up
                echo "Done part 3B"
            fi
        fi

        # PART 4 branch: zero-candidate path still computes the (d-independent) denominator
        if [[ "$PART" == "4" ]]; then
            echo "========== PART 4: Callable genome =========="
            if dnm_candidates_are_empty; then
                echo "[4] 0 LR-validated DNM candidates; d=0. Computing callable genome anyway,"
                echo "[4] then reporting mutation rate = 0."
                write_empty_final_dnmc
            fi
            # Denominator path is independent of the DNM count, so it runs either way.
            calculate_callable_genome
            REfilter_dnm_candidates
            REvalidate_dnmc_with_long_reads
            final_summary
        fi

        

        #############################################
        ## PART CHAIN chain everything from 2b -> 3 -> 3b -> 4 in one go, with dependencies

        if [[ "$PART" == "chain" ]]; then
            echo "========== CHAIN: submit 2b -> 3 -> 3b -> 4 (afterok) =========="
            submit_chain_2b_4_jobs
        fi


        ############################################
        # PART 5 — Fix the blocks that has mismatches
        ############################################
       
        if [[ "$PART" == "5" ]]; then
            echo "========== PART 5: Extract high-mismatch blocks and locally rephase =========="
            rephase_blocks_with_mismatches
        fi
            

        ###############################################
        # PART 6A — Re-run DNM detection with rephased blocks,
        #           compute an INDEPENDENT callable genome for set B,
        #           then build merged file + candidate-weighted mutation rate
        ###############################################

        if [[ "$PART" == "6a" ]]; then
            echo "========== PART 6A: Submit remerge (8h Slurm job) + dependent rest =========="
            submit_rephase_6a_jobs
        fi

        ###############################################
        # PART 6a-merge — INTERNAL: runs inside the merge Slurm job only
        ###############################################
        if [[ "$PART" == "6a-merge" ]]; then
            echo "========== PART 6A-MERGE: remerge all rephased regions =========="
            run_rephase_merge
        fi

        ###############################################
        # PART 6a-rest — INTERNAL: runs inside the dependent Slurm job only
        #   (extract -> count -> filter -> validate -> 6b -> 6d)
        ###############################################
        if [[ "$PART" == "6a-rest" ]]; then
            echo "========== PART 6A-REST: extract + count + filter + validate =========="
            extract_rephase_phased_snp
            count_rephase_mismatches
            filter_rephase_dnm_candidates
            validate_rephase_dnmc_with_parent_bam
            validate_rephase_dnmc_with_hifi_reads

            # Did any rephase candidate survive LR validation?
            REPHASE_B_TSV="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc${REPHASE_SUFFIX}.tsv"
            REPHASE_B_COUNT=0
            if [[ -f "${REPHASE_B_TSV}" ]]; then
                REPHASE_B_COUNT=$(awk 'NR>1 && $1 !~ /^#/ && NF>0' "${REPHASE_B_TSV}" | wc -l)
            fi

            if [[ "${REPHASE_B_COUNT}" -eq 0 ]]; then
                echo "========== PART 6B/6D: No LR-validated rephase DNMs (set B = 0) =========="
                echo "[6a-rest] Skipping callable-genome-rephase; set B contributes 0 DNMs."
                echo "[6a-rest] Final rate = original set-A rate (unchanged)."
                final_summary_rephase_zero
            else
                echo "========== PART 6B: Independent callable genome (rephased blocks) =========="
                calculate_callable_genome_rephase
                REfilter_dnm_candidates_rephase
                REvalidate_dnmc_with_long_reads_rephase

                echo "========== PART 6D: Merge previous + rephase, candidate-weighted rate =========="
                final_summary_with_rephase
            fi
        fi


    done

    echo "Pipeline finished."
}


main "$@"