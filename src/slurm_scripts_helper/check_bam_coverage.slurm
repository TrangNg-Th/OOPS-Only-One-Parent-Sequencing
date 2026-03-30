#!/bin/bash
# Must run in folder of the bam files

#SBATCH -J avg-depth 
#SBATCH -o slm/avg-depth-%j.out
#SBATCH -e slm/avg-depth-%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=nguyetrt@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=8G
#SBATCH -A r00379


## ==============================================
## NOTE
## This script is to check the depth coverage of a bam file using samtools coverage
## It will output the coverage stats for all bam files in the current directory to a file named
## bam.coverage.stat.txt
## ==============================================

set -euo pipefail


# Load modules
module load samtools
echo "Coverage stats for bam files in the current directory:" > bam.coverage.stat.txt

echo "[$(date)] Starting bam coverage check..."
# Check all the available bams in the current directory
for bam in *.bam; do
    echo "Processing $bam file..."
    echo "Bam file: $bam" >> bam.coverage.stat.txt
    samtools coverage $bam >> bam.coverage.stat.txt
    echo "----------------------------------------" >> bam.coverage.stat.txt
done

echo "[$(date)] Done!"