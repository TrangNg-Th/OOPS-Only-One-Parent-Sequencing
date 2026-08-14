#!/bin/bash

#SBATCH -J downsample
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8            # cores given to EACH array task (samtools -@)
#SBATCH --mem=8G
#SBATCH --time=06:00:00             # per BAM; downsampling needs maybe ~4 hrs
#SBATCH --account=r00379            ## Update the account if needed
##SBATCH --partition=CHANGE_ME_IF_YOUR_CLUSTER_REQUIRES_ONE
#SBATCH --output=slurm_output/downspl_%A_%a.out   # %A = array job id, %a = task id
#SBATCH --error=slurm_output/downspl_%A_%a.err
## ---- ONE WORKER PER BAM -------------------------------------------------
## The array runs one independent task per BAM, in parallel.
## Set the upper index to (number of BAMs - 1). With 5 BAMs -> 0-4.
## The optional "%5" caps how many run at once (raise/lower to fit your queue).
#SBATCH --array=0-4%5

set -euo pipefail

mkdir -p slurm_output

## ====================== EDITABLE PARTS ========================
PRJ_DIR="/N/scratch/nguyetrt/OOPS_outputs/original_data/PacBio"

## Bam list -- ONE entry per array index.
## IMPORTANT: keep --array=0-(N-1) above in sync with the length of this list.
LONGREAD_BAMS=(
    "NA12882.CHM13.haplotagged.bam"
    "200082.CHM13.haplotagged.bam"
    "200084.CHM13.haplotagged.bam"
    "NA12877.CHM13.haplotagged.bam"
    "NA12878.CHM13.haplotagged.bam"
)

## FRACTION TO RETAIN
## Change the percentage of sampling by changing the value.
## For example, 10 -> 1.10 means 10% of the reads are kept,
## and 50 -> 1.50 means 50% of the reads are kept.
SAMPLE_FRACTIONS=(27)  # 5 = 5%

# Threads for samtools -- should match --cpus-per-task above
THREADS=8
# =============================================================

# Derived directories -------------------------------------------------
DOWNSAMPLE_DIR="${PRJ_DIR}/downsampled_bams"
mkdir -p "${DOWNSAMPLE_DIR}"

## Load modules
## =============================================================
if command -v module &> /dev/null; then
    module load samtools || true
fi

# --- Pick THIS worker's BAM from the array index --------------------
IDX="${SLURM_ARRAY_TASK_ID}"

if (( IDX >= ${#LONGREAD_BAMS[@]} )); then
    echo "[$(date)] ERROR: array index ${IDX} >= number of BAMs (${#LONGREAD_BAMS[@]}). Fix --array range." >&2
    exit 1
fi

INPUT_BAM_NAME="${LONGREAD_BAMS[$IDX]}"
INPUT_BAM="${PRJ_DIR}/${INPUT_BAM_NAME}"
SAMPLE="${INPUT_BAM_NAME%%.*}"      # e.g. NA12882.CHM13.haplotagged.bam -> NA12882

echo "[$(date)] Task ${IDX} (${SLURM_ARRAY_TASK_ID}) on $(hostname): downsampling ${INPUT_BAM_NAME}"

if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "[$(date)] ERROR: input BAM not found: ${INPUT_BAM}" >&2
    exit 1
fi

for PCT in "${SAMPLE_FRACTIONS[@]}"; do
    # Convert integer percent to samtools format:
    # 50 -> 1.50   25 -> 1.25   10 -> 1.10   5 -> 1.05
    SUBSAMPLE_ARG=$(printf "1.%02d" "${PCT}")

    OUTPUT_BAM="${DOWNSAMPLE_DIR}/${SAMPLE}.CHM13.haplotagged.downsampled_${PCT}pc.bam"

    echo "[$(date)] ${SAMPLE}: downsampling to ${PCT}% -> ${OUTPUT_BAM}"

    samtools view -@ "${THREADS}" -b -s "${SUBSAMPLE_ARG}" "${INPUT_BAM}" -o "${OUTPUT_BAM}"
    samtools index -@ "${THREADS}" "${OUTPUT_BAM}"

    echo "[$(date)] ${SAMPLE}: finished ${PCT}%"
done

echo "[$(date)] Task ${IDX} (${SAMPLE}) completed."
