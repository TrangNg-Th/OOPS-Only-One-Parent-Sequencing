#!/bin/bash

set -eo pipefail
export PS1=${PS1:-}  


################################################################################
# Project: OOPS - Only One Parent Sequencing Mutation Call Pipeline
################################################################################


usage() {
cat << EOF

Only One Parent Sequencing Mutation Call Pipeline
=================================================

This pipeline performs:
  1. HiFi download + Illumina preprocessing + BAM preprocessing
  2. Phasing with long reads + VCF (Whatshap) then haplotagging BAM
  3. De novo mutation (DNM) detection
  4. Local rephasing refinement
  5. Callable genome estimation
  6. Final mutation rate calculation
  7. Optional cleanup of intermediate files

------------------------------------------------------------------
Required arguments:
------------------------------------------------------------------
  --prj-dir           Project root directory
  --sample-child      Child sample ID (e.g. NA12879)
  --sample-parent     Parent sample ID (e.g. NA12878)

------------------------------------------------------------------
Execution control:
------------------------------------------------------------------
  --part <N ...>      Run specific pipeline parts
                      Example: --part 1 2 2b

  --all               Run full pipeline:
                      Parts: 0 1 2 2b 3 3b 4 5

Pipeline Parts:
  0   Just print out the configuration
  1   Data preparation (download + VCF preprocessing + BAM preprocessing)
  2   Initial phasing (Whatshap) + haplotagging cleaned BAM
  2b  First DNM detection
  3   Local rephasing around DNMs
  3b  Refined DNM detection
  4   Callable genome calculation + mutation rate
  5   Remove intermediate files

------------------------------------------------------------------
Optional parameters (with defaults):
------------------------------------------------------------------
  --cpus              CPUs for phasing job (default: 8)
  --reference         Reference file (default: prj-dir/reference/chm13v2.0_maskedY_rCRS.fa)
  --time              Slurm walltime (default: 10:00:00)
  --min-rdepth        Minimum read depth (default: 15)
  --max-rdepth        Maximum read depth (default: 50)
  --gt-qual           Minimum genotype quality (default: 30)
  --nv-quantile       Variant count quantile threshold (default: 0.5)
  --mm-diff-min       Minimum mismatch difference threshold (default: 0.1)
  --min-base-qual     Minimum base quality for LR validation (default: 20)
  --min-map-qual      Minimum mapping quality for LR validation (default: 20)
  --window            Window size around DNM in bp (default: 20000)
  --alt-read-count    Minimum ALT-supporting long reads (default: 8)
  --verbose           T/F verbose long-read validation (default: T)
  --source            Data source configuration file (Example: source.txt  Must be under: ./data/*.txt. default: ./data/source.txt)

------------------------------------------------------------------
Examples:
------------------------------------------------------------------

Run entire pipeline:
  ./pipeline.sh \
    --all \
    --prj-dir /path/to/project \
    --reference /path/tp/project/chm13v2.0_maskedY_rCRS.fa \
    --sample-child NA12879 \
    --sample-parent NA12878

Run only callable genome + final summary:
  ./pipeline.sh \
    --part 4 \
    --prj-dir /path/to/project \
    --reference /path/tp/project/chm13v2.0_maskedY_rCRS.fa \
    --sample-child NA12879 \
    --sample-parent NA12878 \
    --source custom_source.txt

Run data prep + phasing only:
  ./pipeline.sh \
    --part 1 2 \
    --prj-dir /path/to/project \
    --reference /path/tp/project/chm13v2.0_maskedY_rCRS.fa \
    --sample-child NA12879 \
    --sample-parent NA12878

------------------------------------------------------------------

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
            PARTS=("0" "1" "2" "2b" "3" "3b" "4" "5")
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
    SOURCE_FILE="${WORKING_DIR}/data/source.txt"
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

##############################################
# Print configuration summary
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
echo "Alt read count (LR)                   : ${ALT_READ_COUNT}"
echo "Verbose LR validation                 : ${VERBOSE}"
echo "Data source                           : ${SOURCE_FILE}"
echo "Bam file of the child                 : ${BAM_DIR}/${NAME_BAM_CHILD}"
echo "Vcf file                              : ${NAME_VCF}"
echo "Reference file                        : ${REF}"
echo "Parts to run                          : ${PARTS[*]}"
echo "END SUMMARY"
echo "========================================================="
echo ""





##############################################
# Load modules
##############################################
# If on HPS, load these modules, if not, skip
load_modules() {

    if command -v module &> /dev/null; then
        module load bcftools || true
        module load samtools || true
        module load aws-cli || true
        module load conda || true
    fi

}

##############################################
# 1. Download data from AWS
##############################################

download_hifi_data_job() {

    local OUT=./download_hifi_${SAMPLE_CHILD}.slurm
    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J download_hifi_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH -A r00379

set -euo pipefail

module load aws-cli/2.25.5
module load samtools

SAMPLE=${SAMPLE_CHILD}
PRJ_DIR=${PRJ_DIR}
BAM_DIR=${BAM_DIR}
ILLUM_DIR=${ILLUM_DIR}


mkdir -p \${BAM_DIR}
mkdir -p \${ILLUM_DIR}


echo "[download] Downloading HiFi BAM for \${SAMPLE}"

aws s3 cp --no-sign-request ${BAM_CHILD} \${BAM_DIR}/${NAME_BAM_CHILD}

aws s3 cp --no-sign-request ${BAMIDX_CHILD} \${BAM_DIR}/${NAME_BAMIDX_CHILD}

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
    # sbatch "${OUT}"
    rm "${OUT}"
}


##############################################
# 1b. BAM preprocessing
##############################################
generate_bam_preprocessing_job() {

    local OUT=./preprocess_bam_${SAMPLE_CHILD}.slurm

    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J preprocess_bam_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH -A r00379

set -euo pipefail
module load samtools

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

##############################################
# 2. Preprocess VCFs (Extract parent and child vcf file from the full vcf)
##############################################

generate_preprocessing_job() {

    # local DEPENDENCY_JOBID=$1
    local OUT=./preprocess_${SAMPLE_CHILD}.slurm
    # echo "Depending on" ${DEPENDENCY_JOBID}
    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J preprocess_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH -A r00379
#SBATCH -o slurm_output/preprocess_.%j.txt
#SBATCH -e slurm_output/preprocess_.%j_.err


set -euo pipefail
export PS1=\${PS1:-}  

module load bcftools

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
    sbatch "${OUT}"
    rm "${OUT}"
}



###################################################
## 3.5 (optional) Preprocess vcfs ( normalize and compare illumina and hifi vcfs)
###################################################
## These output files are for analysis of False positive. 
## They don't have to be used for phasing as Whatshap will realign things
norm_and_compare_vcfs(){

    local ILLUM_VCF_CHILD_in=${ILLUM_DIR}/${SAMPLE_CHILD}.CHM13.illumina.unphased.noPS.vcf.gz
    local ILLUM_VCF_CHILD_out=${ILLUM_DIR}/${SAMPLE_CHILD}.CHM13.illumina.unphased.noPS.normed.vcf.gz

    # local HIFI_VCF_CHILD_in=${BAM_DIR}/${SAMPLE_CHILD}.CHM13.deepvariant.glnexus.oa.vcf.gz
    # local HIFI_VCF_CHILD_out=${BAM_DIR}/${SAMPLE_CHILD}.CHM13.deepvariant.glnexus.oa.normed.vcf.gz

    local CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

    echo "Normalizing Illumina VCF (autosomes + X,Y)"
    bcftools norm \
      -f ${REF} \
      -m -any \
      -r ${CHRS} \
      ${ILLUM_VCF_CHILD_in} \
      -Oz -o ${ILLUM_VCF_CHILD_out}

    # echo "Normalizing HiFi VCF (autosomes + X,Y)"
    # bcftools norm \
    #   -f ${REF} \
    #   -m -any \
    #   -r ${CHRS} \
    #   ${HIFI_VCF_CHILD_in} \
    #   -Oz -o ${HIFI_VCF_CHILD_out}

    bcftools index ${ILLUM_VCF_CHILD_out}
    # bcftools index ${HIFI_VCF_CHILD_out}
}


##############################################
# 4. Phase the child's illumina vcf using long reads (Slurm job)
##############################################


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
#SBATCH -A r00379
#SBATCH -o slurm_output/build_hapl_%j.txt
#SBATCH -e slurm_output/build_hapl_%j.err

set -euo pipefail
export PS1=\${PS1:-}

module load bcftools
module load samtools
module load conda
conda activate whatshap-env

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

echo "[phase] Checking phased VCF is readable"
bcftools view -h "\${PHASED_VCF}" >/dev/null

echo "[phase] Rebuilding phased VCF index"
rm -f "\${PHASED_VCF}.csi" "\${PHASED_VCF}.tbi"
bcftools index -f "\${PHASED_VCF}"

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

echo "[phase] Done"
EOF

    chmod +x "${OUT}"
    sbatch "${OUT}"
}
##############################################
# 5. Analyze phased mom/child variants output files
##############################################


merge_unphased-parent_phased-child_vcfs() {
    # clean up
    rm -f ./slurm*out
    rm -f ./typescript

    # Make and move phased vcfs to correct directories
    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged

    mkdir -p ${PHASED_VCF}
    mkdir -p ${MERGED_PHASED_VCF}

    
    # Index samples before merging
    bcftools index -f ${ILLUM_DIR}/${SAMPLE_PARENT}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz
    bcftools index -f ${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_LRbam.phased.vcf.gz

    # Merge samples into one vcf
    echo "merging samples"
    bcftools merge -m none -Oz \
        -o ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz \
        ${ILLUM_DIR}/${SAMPLE_PARENT}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz \
        ${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_LRbam.phased.vcf.gz

    bcftools index ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz #index
    echo "done"
}



extract_phased_snp() {
    module load conda

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
        -f '%CHROM,%POS,[%PS,][%GT,][%DP,][%GQ,]\n' \
        > ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_ps.tsv


    ## Fix the extract file to tab delimited format
    conda activate whatshap-env

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


## I think I should apply loose threshold here so that I can get more candidates?
count_shared_alleles_per_PS_block() {
  local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
  local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks
  local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis

  mkdir -p ${MISMATCH_ANALYSIS}
  
  echo " Counting alleles per haplotype block"
  conda activate whatshap-env
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

filter_dnm_candidates() {
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
    FMT/GT[0]==\"0/0\" &&

    FMT/DP>=${MIN_RDEPTH} && FMT/DP<=${MAX_RDEPTH} &&
    FMT/GQ>=${GT_QUAL} &&

    FMT/AD[1:1]>5 &&
    FMT/AD[1:0]>5 &&
    FMT/AD[0:1]==0
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


validate_dnmc_with_hifi_reads(){

  conda activate whatshap-env
  
  local HIFI_BAM="${HP_BAM_PATH}"
  local DNM_BED="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.bed"

  python "${WORKING_DIR}/src/dnmc_readcheck.py" \
    "${SAMPLE_CHILD}" \
    "${HIFI_BAM}" \
    "${DNM_BED}" \
    "${MIN_BASE_QUAL}" \
    "${MIN_MAP_QUAL}" \
    "${WINDW}" \
    "${ALT_READ_COUNT}" \
    "T" \
    "dnmc"
}



### ===================================================================================
## Second validation step
## ======================================================================================

## ======================================================================



##############################################
regenerate_phasing_job() {
    local OUT="${PRJ_DIR}/build_hapl_${SAMPLE_CHILD}.slurm"

    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J build_hapl_${SAMPLE_CHILD}_localphase
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME}
#SBATCH --mem=8G
#SBATCH -A r00379
#SBATCH -o slurm_output/build_hapl_local_%j.txt
#SBATCH -e slurm_output/build_hapl_local_%j.err

set -euo pipefail
export PS1=\${PS1:-}

module load bcftools
module load samtools
module load conda
conda activate whatshap-env

REF="${REF}"
SAMPLE="${SAMPLE_CHILD}"
v="${v}"

ILLUM_VCF="${ILLUM_DIR}/\${SAMPLE}.${NAME_REFERENCE}.illumina.unphased.noPS.vcf.gz"
CLEAN_BAM="${CLEAN_BAM_PATH}"

REGIONS_BED="\${SAMPLE}_dnmc_plusminus\${v}kb.bed"
LOCAL_VCF="\${SAMPLE}_dnmc_plusminus\${v}kb.vcf.gz"
LOCAL_DIPLOID_VCF="\${SAMPLE}_dnmc_plusminus\${v}kb.diploid.vcf.gz"
LOCAL_PHASED_VCF="\${SAMPLE}.illumVCF_LRbam.phased.\${v}kb.vcf.gz"
LOCAL_STATS_TSV="\${SAMPLE}.illumVCF_LRbam.stats.\${v}kb.tsv"
LOCAL_BLOCKS_TSV="\${SAMPLE}.illumVCF_LRbam.blocks.\${v}kb.tsv"

mkdir -p slurm_output

echo "[local phase] Input VCF          : \${ILLUM_VCF}"
echo "[local phase] Clean BAM          : \${CLEAN_BAM}"
echo "[local phase] Regions BED        : \${REGIONS_BED}"
echo "[local phase] Local VCF          : \${LOCAL_VCF}"
echo "[local phase] Local diploid VCF  : \${LOCAL_DIPLOID_VCF}"
echo "[local phase] Local phased VCF   : \${LOCAL_PHASED_VCF}"

if [[ ! -f "\${ILLUM_VCF}" ]]; then
    echo "ERROR: Illumina VCF not found: \${ILLUM_VCF}"
    exit 1
fi

if [[ ! -f "\${CLEAN_BAM}" ]]; then
    echo "ERROR: Clean BAM not found: \${CLEAN_BAM}"
    exit 1
fi

if [[ ! -f "\${REGIONS_BED}" ]]; then
    echo "ERROR: Regions BED not found: \${REGIONS_BED}"
    exit 1
fi

if [[ ! -f "\${CLEAN_BAM}.bai" ]]; then
    echo "[local phase] BAM index not found. Creating index for clean BAM."
    samtools index "\${CLEAN_BAM}"
fi

echo "[local phase] Removing stale local outputs"
rm -f "\${LOCAL_VCF}" "\${LOCAL_VCF}.csi" "\${LOCAL_VCF}.tbi"
rm -f "\${LOCAL_DIPLOID_VCF}" "\${LOCAL_DIPLOID_VCF}.csi" "\${LOCAL_DIPLOID_VCF}.tbi"
rm -f "\${LOCAL_PHASED_VCF}" "\${LOCAL_PHASED_VCF}.csi" "\${LOCAL_PHASED_VCF}.tbi"
rm -f "\${LOCAL_STATS_TSV}" "\${LOCAL_BLOCKS_TSV}"

echo "[local phase] Subsetting Illumina VCF to local regions"
bcftools view \
  -R "\${REGIONS_BED}" \
  -s "\${SAMPLE}" \
  -Oz -o "\${LOCAL_VCF}" \
  "\${ILLUM_VCF}"

echo "[local phase] Indexing local VCF"
bcftools index -f "\${LOCAL_VCF}"

N_LOCAL=\$(bcftools index -n "\${LOCAL_VCF}")
echo "[local phase] Number of records in local VCF: \${N_LOCAL}"

if [[ "\${N_LOCAL}" -eq 0 ]]; then
    echo "ERROR: No variants found in local regions: \${LOCAL_VCF}"
    exit 1
fi

echo "[local phase] Building diploid-only local VCF to avoid mixed-ploidy crashes"
bcftools view -Ou "\${LOCAL_VCF}" \
  | bcftools view -i 'GT="het" || GT="hom" || GT="ref" || GT="alt"' \
  -Oz -o "\${LOCAL_DIPLOID_VCF}"

echo "[local phase] Indexing diploid-only local VCF"
bcftools index -f "\${LOCAL_DIPLOID_VCF}"

N_DIPLOID=\$(bcftools index -n "\${LOCAL_DIPLOID_VCF}")
echo "[local phase] Number of diploid local records retained: \${N_DIPLOID}"

if [[ "\${N_DIPLOID}" -eq 0 ]]; then
    echo "ERROR: No diploid variants remained after filtering: \${LOCAL_DIPLOID_VCF}"
    exit 1
fi

echo "[local phase] Optional sanity check: count records whose GT is not diploid-style"
bcftools query -f '[%GT\n]' "\${LOCAL_DIPLOID_VCF}" | awk '
    \$1 !~ /^[.0-9]+[\/|][.0-9]+$/ {bad++}
    END {print "[local phase] Non-diploid-style GT records remaining:", bad+0}
'

echo "[local phase] Running whatshap on diploid-only local VCF"
whatshap phase \
  --reference "\${REF}" \
  --ignore-read-groups \
  -o "\${LOCAL_PHASED_VCF}" \
  "\${LOCAL_DIPLOID_VCF}" "\${CLEAN_BAM}"

echo "[local phase] Checking phased local VCF is readable"
bcftools view -h "\${LOCAL_PHASED_VCF}" >/dev/null

echo "[local phase] Rebuilding phased local VCF index"
rm -f "\${LOCAL_PHASED_VCF}.csi" "\${LOCAL_PHASED_VCF}.tbi"
bcftools index -f "\${LOCAL_PHASED_VCF}"

echo "--ignore-read-groups"
whatshap stats \
  --tsv="\${LOCAL_STATS_TSV}" \
  --block-list="\${LOCAL_BLOCKS_TSV}" \
  "\${LOCAL_PHASED_VCF}"

echo "[local phase] Done"
EOF

    chmod +x "${OUT}"
    sbatch "${OUT}"
}



remerge_unphased-parent_phased-child_vcfs() {
    echo "Re-merging parent vcf and child vcf"
    # clean up
    rm -f ./slurm*out
    rm -f ./typescript

    # Make and move phased vcfs to correct directories
    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged


    mkdir -p ${PHASED_VCF}
    mkdir -p ${MERGED_PHASED_VCF}

    # move child vcf to correct place
    ls ${WORKING_DIR}/${SAMPLE_CHILD}.illumVCF_LRbam* 1>/dev/null 2>&1 && \
        mv -f ${WORKING_DIR}/${SAMPLE_CHILD}.illumVCF_LRbam* ${PHASED_VCF}/
    
    # Index samples before merging
    bcftools index -f ${ILLUM_DIR}/${SAMPLE_PARENT}.CHM13.illumina.unphased.noPS.vcf.gz
    bcftools index -f ${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_LRbam.phased.${v}kb.vcf.gz

    # Merge samples into one vcf
    bcftools merge -m none -Oz \
        -o ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.${v}kb.vcf.gz \
        ${ILLUM_DIR}/${SAMPLE_PARENT}.CHM13.illumina.unphased.noPS.vcf.gz \
        ${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_LRbam.phased.${v}kb.vcf.gz

    bcftools index ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.${v}kb.vcf.gz #index
    echo "done"
}


reextract_phased_snp() {
    local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged
    local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks

    module load bcftools 
    module load conda

    echo "Reextracting phased snps"
    

    mkdir -p ${EXTRACT_HAPLBLOCK}
    mkdir -p ${MERGED_PHASED_VCF}
    mkdir -p ${EXTRACT_HAPLBLOCK}

    bcftools view ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.${v}kb.vcf.gz \
    -i 'TYPE="snp" && N_ALT=1 && GT[0]!="mis" && GT[1]!="mis" && FMT/PS[1]!="."' \
    -Ou \
    | bcftools query \
        -s ${SAMPLE_PARENT},${SAMPLE_CHILD} \
        -f '%CHROM,%POS,[%PS,][%GT,][%DP,][%GQ,]\n' \
        > ${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps.tsv

 
    ## fix again the extracted format into tab delim file
    conda activate whatshap-env
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
  
  conda activate whatshap-env
  python "${WORKING_DIR}/src/count_mismatches.py" \
  "${EXTRACT_HAPLBLOCK}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_${v}kb.ps.tsv" \
  "${MISMATCH_ANALYSIS}" \
  "${MIN_RDEPTH}" \
  "${MAX_RDEPTH}" \
  "${SAMPLE_PARENT}" \
  "${SAMPLE_CHILD}" \
  "${GT_QUAL}" \
    0 \
    0 \
    "${WINDW}" \
    "${RECOUNT}" \
    "${DNMC_FILE}"

    echo "See outputs in ${MISMATCH_ANALYSIS}"
    echo "Done"

}

clean_up(){
  local SLICE_DIR=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis/sliced/
  
  mkdir -p ${SLICE_DIR}

  mv ${SAMPLE_CHILD}*plusminus* ${SLICE_DIR}
  mv ${SAMPLE_CHILD}*LR_validated* ${SLICE_DIR}

}


calculate_callable_genome() {
  local PHASED_VCF=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
  local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks
  local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis
  local POST_ANALYSIS=${PHASED_VCF}/mismatch_analysis/denum_calcul
  
  mkdir -p ${POST_ANALYSIS}
  conda activate whatshap-env
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


REvalidate_dnmc_with_hifi_reads(){

  conda activate whatshap-env
  
  local HIFI_BAM="${HP_BAM_PATH}"
  local DNM_BED="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis/denum_calcul/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.bed"
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
    "hetc"
}

final_summary() {

    local WORKING_DIR="${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf"
    local DENUM_DIR="${WORKING_DIR}/mismatch_analysis/denum_calcul"
    local dnmc_file="${WORKING_DIR}/final_dnmc_${SAMPLE_CHILD}-from-${SAMPLE_PARENT}.tsv"
    local nb_qualified_snps="${DENUM_DIR}/${SAMPLE_CHILD}_LR_validated_hetc.bed"
    local callable_genome_file="${DENUM_DIR}/callable_genome.txt"

    mkdir -p "${DENUM_DIR}"

    compgen -G "./${SAMPLE_CHILD}*plusminus*" > /dev/null && \
        mv ./${SAMPLE_CHILD}*plusminus* "${DENUM_DIR}/"

    compgen -G "./${SAMPLE_CHILD}*LR_validated*" > /dev/null && \
        mv ./${SAMPLE_CHILD}*LR_validated* "${DENUM_DIR}/"

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

final_cleanup(){
    
    local WORKING_DIR2=${PRJ_DIR}/${SAMPLE_CHILD}_phasedvcf
    
    # Remove all intermediate files
    rm ${WORKING_DIR2}/${SAMPLE_CHILD}.illumVCF_LRbam*
    rm ${WORKING_DIR2}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_mismatch*
    rm ${WORKING_DIR2}/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc*

}

##################################################
# Main
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
        # PART 1 — Data Preparation
        ############################################
        if [[ "$PART" == "1" ]]; then
            echo "========== PART 1: Data preparation =========="

            # Download data from AWS
            # DOWNLOAD_JOBID=$(download_hifi_data_job)
            # echo "Download job submitted: ${DOWNLOAD_JOBID}"

     
            PREPROCESS_VCF_JOBID=$(generate_preprocessing_job)
            echo "VCF preprocessing job submitted: ${PREPROCESS_VCF_JOBID}"

            PREPROCESS_BAM_JOBID=$(generate_bam_preprocessing_job)
            echo "BAM preprocessing (check HP tag) job submitted: ${PREPROCESS_BAM_JOBID}" 


            # Normalize two sources of vcf in the child ( Illumina / Hifi) - optional
            # norm_and_compare_vcfs
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
     

            # Additional filters with hifi suppport
            echo "Validating DNM candidates with HiFi read support"
            validate_dnmc_with_hifi_reads

        fi


        ############################################
        # PART 3 — Local rephasing
        ############################################
        if [[ "$PART" == "3" ]]; then
            echo "========== PART 3: Local rephasing =========="
            # Rephase the areas around the dnmc
            regenerate_phasing_job
    
        fi


        ############################################
        # PART 3B — Re-merge + refined DNM
        ############################################
        if [[ "$PART" == "3b" ]]; then
            echo "========== PART 3B: Refined DNM detection =========="

            # remerge the phased dnmc with the parent's vcf
            remerge_unphased-parent_phased-child_vcfs

            # Reextract the phased dnmc
            reextract_phased_snp

            # Recount within a block ==> Final results
            recount_shared_alleles_per_PS_block

            # Cleanup
            clean_up

            echo "Done part 3B"
        fi


        ############################################
        # PART 4 — Callable genome
        ############################################
        if [[ "$PART" == "4" ]]; then
            echo "========== PART 4: Callable genome =========="
            
            # Find the denominator of the mutation rate
            calculate_callable_genome

            # Use the filters applied to the dnm to apply to the randomly selected snps
            REfilter_dnm_candidates

            # Apply another filter - long read filter to those snps ==> see how many snps can we obtain
            REvalidate_dnmc_with_hifi_reads

            # Give a final summary results
            final_summary
        fi


        ############################################
        # PART 5 — Clean up everything
        ############################################
        if [[ "$PART" == "5" ]]; then
            echo "========== PART 5: Removing all the intermediate files  =========="
            
           # Clean up all intermediate files
           final_cleanup
        fi


    done

    echo "Pipeline finished."
}


main "$@"