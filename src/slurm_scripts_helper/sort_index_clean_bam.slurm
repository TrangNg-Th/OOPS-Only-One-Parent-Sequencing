#!/bin/bash
# Must run from the repo root (so relative paths / slurm_output exist), OR
# just point it straight at a *_clean.bam with a positional argument.

#SBATCH -J sortidx-clean
#SBATCH -o slurm_output/sortidx-clean-%j.txt
#SBATCH -e slurm_output/sortidx-clean-%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=CHANGE_ME_TO_YOUR_EMAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --account=r00379 


set -euo pipefail

## =============================================================
## PURPOSE
## Finish the tail of PART 1B (BAM preprocessing) for ONE clean BAM when the
## original slurm job ran out of time before "sort clean BAM" / "index clean
## BAM" (main.sh lines ~560-565). It (re)sorts the *_clean.bam by coordinate
## and (re)builds its .bai. HP-tag stripping is assumed already done — this
## only redoes the sort + index, so it is safe to re-run.
##
## USAGE
##   # A) Point it directly at the clean BAM:
##   sbatch src/slurm_scripts_helper/sort_index_clean_bam.slurm \
##          /N/scratch/nguyetrt/OOPS_outputs/OOPS_Pacbio_10x_78_91/bam/<child>_clean.bam
##
##   # B) Or let it derive the path the same way main.sh does, by editing
##   #    PRJ_DIR + SOURCE_FILE below and running with no argument:
##   sbatch src/slurm_scripts_helper/sort_index_clean_bam.slurm
##
##   # Runs fine with plain bash too (no slurm):
##   bash src/slurm_scripts_helper/sort_index_clean_bam.slurm <clean_bam>
## =============================================================

## ---- Edit these ONLY if you are NOT passing the BAM path as an argument ----
PRJ_DIR="/N/scratch/nguyetrt/OOPS_outputs/OOPS_Pacbio_10x_78_91"
SOURCE_FILE="/N/project/mutation_rate_Mmulatta/CEPH_data/OOPS-Only-One-Parent-Sequencing/data/source_Pacbio_10x_78_91.txt"
## ---------------------------------------------------------------------------

## Load modules
if command -v module &> /dev/null; then
    module load samtools || true
fi

THREADS="${SLURM_CPUS_PER_TASK:-4}"

## Resolve the clean BAM path
## -------------------------------------------------------------
if [[ -n "${1:-}" ]]; then
    CLEAN_BAM="$1"
else
    # Mirror main.sh: BAM_DIR=${PRJ_DIR}/bam ,
    # CLEAN_BAM=${BAM_DIR}/$(basename NAME_BAM_CHILD .bam)_clean.bam
    if [[ ! -f "${SOURCE_FILE}" ]]; then
        echo "ERROR: source file not found: ${SOURCE_FILE}" >&2
        echo "       Pass the clean BAM path as an argument instead." >&2
        exit 1
    fi
    NAME_BAM_CHILD="$(grep -v '^##' "${SOURCE_FILE}" \
                      | grep -E '^[[:space:]]*NAME_BAM_CHILD=' \
                      | head -n1 | cut -d= -f2- | tr -d '[:space:]')"
    if [[ -z "${NAME_BAM_CHILD}" ]]; then
        echo "ERROR: NAME_BAM_CHILD not found in ${SOURCE_FILE}" >&2
        exit 1
    fi
    BAM_DIR="${PRJ_DIR}/bam"
    BASENAME="$(basename "${NAME_BAM_CHILD}" .bam)"
    CLEAN_BAM="${BAM_DIR}/${BASENAME}_clean.bam"
fi

echo "[$(date)] sort+index target: ${CLEAN_BAM} (threads=${THREADS})"

if [[ ! -f "${CLEAN_BAM}" ]]; then
    echo "ERROR: clean BAM not found: ${CLEAN_BAM}" >&2
    echo "       (expected the *_clean.bam produced by PART 1B)" >&2
    exit 1
fi

## Verify the clean BAM itself is not truncated before we touch it
samtools quickcheck "${CLEAN_BAM}" || {
    echo "ERROR: ${CLEAN_BAM} fails samtools quickcheck (truncated/corrupt)." >&2
    echo "       PART 1B likely died mid-write; re-run PART 1B from the original BAM." >&2
    exit 1
}

## Remove any leftover temp from an interrupted sort, then sort -> index.
## Sort writes to a temp and is only mv'd into place after quickcheck passes,
## so an interruption never leaves a half-written clean BAM.
rm -f "${CLEAN_BAM}.sorted" "${CLEAN_BAM}".tmp.*.bam 2>/dev/null || true

echo "[$(date)] Sorting clean BAM"
samtools sort -@ "${THREADS}" -o "${CLEAN_BAM}.sorted" "${CLEAN_BAM}"
samtools quickcheck "${CLEAN_BAM}.sorted"
mv "${CLEAN_BAM}.sorted" "${CLEAN_BAM}"

echo "[$(date)] Indexing clean BAM"
samtools index -@ "${THREADS}" "${CLEAN_BAM}"

if [[ ! -f "${CLEAN_BAM}.bai" ]]; then
    echo "ERROR: index not created: ${CLEAN_BAM}.bai" >&2
    exit 1
fi

echo "[$(date)] Done. Sorted + indexed: ${CLEAN_BAM}"
echo "         You can now run PART 2 for this pair."
