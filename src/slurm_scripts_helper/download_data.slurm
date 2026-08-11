#!/bin/bash
#SBATCH --job-name=download_data
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --account=r00379                       # <-- change to your Slurm account
#SBATCH --output=slurm_output/download_data_%j.out
#SBATCH --error=slurm_output/download_data_%j.err
## If your cluster requires a partition/queue, uncomment and set it:
## SBATCH --partition=general
#
# NOTE: #SBATCH directives are read by sbatch BEFORE the shell runs, so they
# cannot contain shell variables. Keep them literal, or override at submission
# time, e.g.:  sbatch --job-name=download_NA12879 --account=myacct download_data.slurm
#
# The log directory in --output/--error must already exist at submission time:
#   mkdir -p slurm_output
# before you run sbatch, or the job will fail to start.

set -euo pipefail

# ============================ EDIT THESE ============================
PRJ_DIR="/N/scratch/nguyetrt/OOPS_outputs/OOPS_Pacbio_10x_79_77"
SAMPLE_CHILD="NA12879"
SAMPLE_PARENT="NA12877"
SOURCE_FILE="/N/project/mutation_rate_Mmulatta/CEPH_data/OOPS-Only-One-Parent-Sequencing/data/source_Pacbio_10x_79_77.txt"
# ===================================================================

# ---- Load the data-source config -----------------------------------------
# Provides: BAM_CHILD_URL, VCF_PATH, BAM_PARENT_URL, REFERENCE_PATH,
#           NAME_BAM_CHILD, NAME_BAMIDX_CHILD, NAME_VCF_FILE,
#           NAME_REFERENCE_FILE, NAME_REFERENCE, NAME_BAM_PARENT,
#           NAME_BAMIDX_PARENT   ('##' comment lines are skipped)
if [[ ! -f "${SOURCE_FILE}" ]]; then
    echo "ERROR: source config not found: ${SOURCE_FILE}" >&2
    exit 1
fi
set -a
source <(grep -v '^##' "${SOURCE_FILE}")
set +a

# ---- Modules (HPC only; harmless elsewhere) ------------------------------
if command -v module &> /dev/null; then
    module load aws-cli/2.25.5 || true
    module load samtools       || true
fi

# ---- Derived directories and paths (mirror main.sh) ----------------------
BAM_DIR="${PRJ_DIR}/bam"
ILLUM_DIR="${PRJ_DIR}/illumina-dragen"
REF_DIR="${PRJ_DIR}/reference"

BAM_CHILD="${BAM_CHILD_URL}/${NAME_BAM_CHILD}"
BAMIDX_CHILD="${BAM_CHILD_URL}/${NAME_BAMIDX_CHILD}"
VCF_FULL_PATH="${VCF_PATH}/${NAME_VCF_FILE}"
NAME_VCF="${ILLUM_DIR}/${NAME_VCF_FILE}"
REF_GZ_NAME="$(basename "${REFERENCE_PATH}")"   # e.g. chm13v2.0_maskedY_rCRS.fa.gz

mkdir -p "${BAM_DIR}" "${ILLUM_DIR}" "${REF_DIR}"

# ---- Child long-read BAM -------------------------------------------------
echo "[download] Child long-read BAM for ${SAMPLE_CHILD}"
aws s3 cp --no-sign-request "${BAM_CHILD}"    "${BAM_DIR}/${NAME_BAM_CHILD}"
aws s3 cp --no-sign-request "${BAMIDX_CHILD}" "${BAM_DIR}/${NAME_BAMIDX_CHILD}"

# ---- Parent Illumina BAM -------------------------------------------------
echo "[download] Parent Illumina BAM for ${SAMPLE_PARENT}"
aws s3 cp --no-sign-request "${BAM_PARENT_URL}/${NAME_BAM_PARENT}"    "${BAM_DIR}/${NAME_BAM_PARENT}"
aws s3 cp --no-sign-request "${BAM_PARENT_URL}/${NAME_BAMIDX_PARENT}" "${BAM_DIR}/${NAME_BAMIDX_PARENT}"

if [[ ! -f "${BAM_DIR}/${NAME_BAM_PARENT}.bai" \
   && ! -f "${BAM_DIR}/$(basename "${NAME_BAM_PARENT}" .bam).bai" ]]; then
    echo "[download] Parent BAM index not found; building one"
    samtools index "${BAM_DIR}/${NAME_BAM_PARENT}"
fi

# ---- Reference (download .fa.gz, decompress, index) ----------------------
echo "[download] Reference file"
if [[ ! -f "${REF_DIR}/${NAME_REFERENCE_FILE}" ]]; then
    wget -P "${REF_DIR}" "${REFERENCE_PATH}"
    gunzip -f "${REF_DIR}/${REF_GZ_NAME}"       # -> ${NAME_REFERENCE_FILE}
fi
samtools faidx "${REF_DIR}/${NAME_REFERENCE_FILE}"

# ---- Illumina joint VCF --------------------------------------------------
echo "[download] Illumina VCF"
aws s3 cp --no-sign-request "${VCF_FULL_PATH}" "${NAME_VCF}"

echo "[download] Done"
