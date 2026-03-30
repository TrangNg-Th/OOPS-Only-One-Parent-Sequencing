#!/bin/bash
#SBATCH -J addRG_NA12879
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH -A r00379
#SBATCH -o logs/addRG_%x_%j.out
#SBATCH -e logs/addRG_%x_%j.err

set -euo pipefail


# =============================================
## This script should be executed inside the same directory as the
## bam files of ONT. 
## The problem is that in Oxford Nanopore bam file : There are some that don't have READ GROUP header
## meaning, there are missing sample names
## This script adds read group to the bam file and then indexes it.
## Note: This is not necessary for the analysis, but it is good practice to have read groups in the bam file.

# =============================================
# Change this for other samples if needed
# =============================================
SAMPLE=NA12879
INPUT_BAM=${SAMPLE}.CHM13.minimap2.bam
OUTPUT_BAM=${SAMPLE}.CHM13.minimap2.fixed.bam

# =============================================
# Load modules
# =============================================
module load samtools

# =============================================
# Add Read Group
# =============================================
echo "[$(date)] Adding read group to ${INPUT_BAM}"

samtools addreplacerg \
    -r "ID:1 SM:${SAMPLE} PL:ONT" \
    -o ${OUTPUT_BAM} \
    ${INPUT_BAM}

# =============================================
# Index BAM
# =============================================
echo "[$(date)] Indexing ${OUTPUT_BAM}"

samtools index ${OUTPUT_BAM}

echo "[$(date)] Done!"