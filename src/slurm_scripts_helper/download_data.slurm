#!/bin/bash
#SBATCH --job-name=download_data
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --account=r00379                      ## Update the account if neede
#SBATCH --output=slurm_output/download_data_%j.out
#SBATCH --error=slurm_output/download_data_%j.err


## Output directories for slurm logs
mkdir -p slurm_output


set -euo pipefail

# ============================ EDITABLE PARTS  ============================
PRJ_DIR="/N/scratch/nguyetrt/OOPS_outputs/original_data"
mkdir -p "${PRJ_DIR}"

# S3 base URLs (one per file type) 
BAM_LONGREAD_URL="s3://platinum-pedigree-data/data/hifi/mapped/CHM13"       # long-read (HiFi/ONT) BAMs
BAM_ILLUMINA_URL="s3://platinum-pedigree-data/data/illumina/mapped/CHM13"   # Illumina BAMs
VCF_PATH="s3://platinum-pedigree-data/variants/small_variants/CHM13"        # joint VCF location

# Reference (single) 
REFERENCE_PATH="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz"
NAME_REFERENCE_FILE="chm13v2.0_maskedY_rCRS.fa"   # decompressed name (what the pipeline uses)

# Joint VCF (single) 
NAME_VCF_FILE="CEPH1463.CHM13.illumina-dragen.oa.vcf.gz"

# BAM lists (FULL object names; one entry per file) -------------------
LONGREAD_BAMS=(
    "NA12879.CHM13.haplotagged.downsampled_27pc.bam"
)

ILLUMINA_BAMS=(
    "NA12877.illumina.bam"
)


# Derived directories -------------------------------------------------
BAM_DIR="${PRJ_DIR}/PacBio"
ILLUM_DIR="${PRJ_DIR}/illumina-dragen"
REF_DIR="${PRJ_DIR}/reference"
REF_GZ_NAME="$(basename "${REFERENCE_PATH}")"     # e.g. chm13v2.0_maskedY_rCRS.fa.gz

mkdir -p "${BAM_DIR}" "${ILLUM_DIR}" "${REF_DIR}"

# ===================================================================

# Modules (HPC only; harmless elsewhere) ------------------------------
if command -v module &> /dev/null; then
    module load aws-cli/2.25.5 || true
    module load samtools       || true
fi



# Helper: download one BAM + its index --------------------------------
download_bam () {

    ## 
    local base="$1" bam="$2"
    local dest="${BAM_DIR}/${bam}"

    echo "[download]   ${bam}"
    aws s3 cp --no-sign-request "${base}/${bam}" "${dest}"


    local bai_a="${bam}.bai"                        # file.bam.bai
    local bai_b="$(basename "${bam}" .bam).bai"     # file.bai

    # Download bam and its index if possible, if not, just use samtoosl
    if aws s3 cp --no-sign-request "${base}/${bai_a}" "${BAM_DIR}/${bai_a}" 2>/dev/null; then
        :
    elif aws s3 cp --no-sign-request "${base}/${bai_b}" "${BAM_DIR}/${bai_b}" 2>/dev/null; then
        :
    elif [[ -f "${BAM_DIR}/${bai_a}" || -f "${BAM_DIR}/${bai_b}" ]]; then
        :   # already present locally
    else
        echo "[download]   index not in bucket; building with samtools index"
        samtools index "${dest}"
    fi
}

# Long-read BAMs ------------------------------------------------------
echo "[download] Long-read BAMs (${#LONGREAD_BAMS[@]})"
for bam in "${LONGREAD_BAMS[@]}"; do
    download_bam "${BAM_LONGREAD_URL}" "${bam}"
    echo "[download]   ${bam} done"
done

# Illumina BAMs -------------------------------------------------------
echo "[download] Illumina BAMs (${#ILLUMINA_BAMS[@]})"
for bam in "${ILLUMINA_BAMS[@]}"; do
    download_bam "${BAM_ILLUMINA_URL}" "${bam}"
    echo "[download]   ${bam} done"
done

# Reference (download .fa.gz, decompress, index) ----------------------
echo "[download] Reference file"
if [[ ! -f "${REF_DIR}/${NAME_REFERENCE_FILE}" ]]; then
    wget -P "${REF_DIR}" "${REFERENCE_PATH}"
    gunzip -f "${REF_DIR}/${REF_GZ_NAME}"       # -> ${NAME_REFERENCE_FILE}
fi
samtools faidx "${REF_DIR}/${NAME_REFERENCE_FILE}"

# Joint VCF -----------------------------------------------------------
echo "[download] Joint VCF"
aws s3 cp --no-sign-request "${VCF_PATH}/${NAME_VCF_FILE}" "${ILLUM_DIR}/${NAME_VCF_FILE}"

echo "[download] Done"
