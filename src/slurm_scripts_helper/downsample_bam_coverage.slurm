#!/bin/bash

#SBATCH -J downspl 
#SBATCH -o slurm_output/1.downspl.%j.txt
#SBATCH -e slurm_output/1.downspl.%j_.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nguyetrt@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00 # downsampling need maybe 4 hrs
#SBATCH --mem=8G
#SBATCH -A r00379


set -euo pipefail

# =============================================================
## This script should be executed inside the same directory as the bam file to be downsampled.
## =============================================================



## Load moduless
## =============================================================
module load samtools


## Change the percentage of sampling by changing the value after -s. 
## For example, -s 1.10 means 10% of the reads will be kept, 
## and -s 1.50 means 50% of the reads will be kept.
## =============================================================
SAMPLE=NA12879
INPUT_BAM=${SAMPLE}.CHM13.haplotagged.bam

# Percentages to retain
SAMPLE_FRACTIONS=(5)  # 5 = 5%

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
