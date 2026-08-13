#!/bin/bash

#SBATCH -J downsample 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o slurm_output/1.downspl.%j.txt
#SBATCH -e slurm_output/1.downspl.%j_.err
#SBATCH --time=23:00:00 # downsampling need maybe 4 hrs
#SBATCH --account=r00379                      ## Update the account if neede
#SBATCH --output=slurm_output/downspl_%j.out
#SBATCH --error=slurm_output/downspl_%j.err


set -euo pipefail

## ====================== EDITABLE PARTS ========================
PRJ_DIR="/N/scratch/nguyetrt/OOPS_outputs/original_data/PacBio"

## Bam list
LONGREAD_BAMS=(
    "NA12882.CHM13.haplotagged.bam"
    "200082.CHM13.haplotagged.bam"
    "200084.CHM13.haplotagged.bam"
    "NA12877.CHM13.haplotagged.bam"
    "NA12878.CHM13.haplotagged.bam"

)

## FRACTION TO RETAIN
## Change the percentage of sampling by changing the value after -s. 
## For example, -s 1.10 means 10% of the reads will be kept, 
## and -s 1.50 means 50% of the reads will be kept.
SAMPLE_FRACTIONS=(5)  # 5 = 5%

# Derived directories -------------------------------------------------
DOWNSAMPLE_DIR="${PRJ_DIR}/downsampled_bams"

mkdir -p "${DOWNSAMPLE_DIR}"

## Load moduless
## =============================================================
if command -v module &> /dev/null; then
    module load samtools || true
fi

# Helper : downsample a single bam file


echo "[$(date)] Starting downsampling for ${INPUT_BAM}"

for PCT in "${SAMPLE_FRACTIONS[@]}"; do
    # Convert integer percent to samtools format:
    # 50 -> 1.50
    # 25 -> 1.25
    # 10 -> 1.10
    SUBSAMPLE_ARG=$(printf "1.%02d" "${PCT}")

    OUTPUT_BAM="${SAMPLE}.CHM13.haplotagged.downsampled_${PCT}pc.bam"

    echo "[$(date)] Downsampling to ${PCT}% -> ${OUTPUT_BAM}"

    samtools view -@ 8 -b -s "${SUBSAMPLE_ARG}" "${INPUT_BAM}" -o "${OUTPUT_BAM}"
    samtools index -@ 8 "${OUTPUT_BAM}"

    echo "[$(date)] Finished ${PCT}%"
done

echo "[$(date)] All downsampling jobs completed."
