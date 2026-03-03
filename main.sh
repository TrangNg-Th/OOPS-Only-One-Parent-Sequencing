#!/bin/bash
export PS1=${PS1:-}
set -euo pipefail


################################################################################
# Project: Platinum pedigree analysis
################################################################################

##############################################
# Configuration (user input)
##############################################
usage() {
  cat << EOF

  Only One Parent Sequencing Mutation Call Pipeline
  =====================================================

  Required arguments : 
    --prj-dir           Project root directory
    --sample-child      Child sample ID in the VCF (e.g. NA12879)
    --sample-parent     Parent sample ID in the VCF (e.g. NA12878)

  Optional arguments (with defaults):
    --cpus              Number of CPUs for phasing job (default: 8)
    --time              Slurm walltime (default: 4:00:00)
    --min-rdepth        Minimum read depth (default: 15)
    --max-rdepth        Maximum read depth (default: 50)
    --gt-qual           Minimum genotype quality (default: 30)
    --nv-quantile       Noise quantile threshold (default: 0.75)
    --mm-diff-min       Minimum mismatch difference threshold (default: 0.1)
    --min-base-qual     Minimum base quality for LR validation (default: 20)
    --min-map-qual      Minimum mapping quality for LR validation (default: 20)
    --window            Window size around DNM (default: 20000)
    --alt-read-count    Minimum ALT supporting long reads (default: 8)
    --verbose           T/F for verbose LR validation (default: T)

  Example:
    ./pipeline.sh \
      --prj-dir /N/project/mutation_rate_Mmulatta/platinum-ped-data/aws-data \
      --sample-child NA12879 \
      --sample-parent NA12878

EOF
  exit 1
}

##############################################
# Default values
##############################################

CPUS=8
TIME="4:00:00"
MIN_RDEPTH=15
MAX_RDEPTH=50
GT_QUAL=30
NV_QUANTILE=0.75
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
            PARTS=("1" "2" "2b" "3" "3b" "4")
            shift
            ;;

        --prj-dir) PRJ_DIR="$2"; shift 2 ;;
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

ILLUM_DIR=${PRJ_DIR}/variants/small_variants/illumina-dragen
HIFI_DIR=${PRJ_DIR}/variants/small_variants/hifi
REF=${PRJ_DIR}/reference/chm13v2.0_maskedY_rCRS.fa
v=$((WINDW / 1000))

##############################################
# Print configuration summary
##############################################

echo "================ Pipeline Configuration ================"
echo "Project directory        : ${PRJ_DIR}"
echo "Child sample             : ${SAMPLE_CHILD}"
echo "Parent sample            : ${SAMPLE_PARENT}"
echo "CPUs                     : ${CPUS}"
echo "Walltime                 : ${TIME}"
echo "Min read depth           : ${MIN_RDEPTH}"
echo "Max read depth           : ${MAX_RDEPTH}"
echo "Genotype quality cutoff  : ${GT_QUAL}"
echo "Noise quantile           : ${NV_QUANTILE}"
echo "Mismatch diff threshold  : ${MM_DIFF_MIN}"
echo "Min base quality (LR)    : ${MIN_BASE_QUAL}"
echo "Min mapping quality (LR) : ${MIN_MAP_QUAL}"
echo "Window size (bp)         : ${WINDW}"
echo "Window size (kb)         : ${v}"
echo "Alt read count (LR)      : ${ALT_READ_COUNT}"
echo "Verbose LR validation    : ${VERBOSE}"
echo "========================================================="
echo ""


## =================================================



PRJ_DIR="/N/project/mutation_rate_Mmulatta/platinum-ped-data/aws-data"
CPUS=8
TIME=4:00:00
SAMPLE_CHILD=NA12879
SAMPLE_PARENT=NA12878 ## 78: mom, 77: dad
MIN_RDEPTH=15
MAX_RDEPTH=50
ILLUM_DIR=${PRJ_DIR}/variants/small_variants/illumina-dragen
REF=${PRJ_DIR}/reference/chm13v2.0_maskedY_rCRS.fa
HIFI_DIR=${PRJ_DIR}/variants/small_variants/hifi
GT_QUAL=30
NV_QUANTILE=0.75
MM_DIFF_MIN=0.1
MIN_BASE_QUAL=20
MIN_MAP_QUAL=20
WINDW=20000
v=$((WINDW / 1000))

RECOUNT=T
NOTRECOUNT=F ## give alternative option to count mismatches 
VERBOSE=T  # debugging print for dnmc_readcheck.py to exploit LR bam
ALT_READ_COUNT=8 # number of alternative read count on long read data


##############################################
# Load modules
##############################################
load_modules() {
    module load aws-cli/2.25.5
    module load bcftools
    module load conda
}



##############################################
# 1. Download data from AWS
##############################################
download_data() {

    # HiFi BAMs
    ## ============================================================================================
    ## CHILD - REF CHM13
    aws s3 cp --no-sign-request \
      s3://platinum-pedigree-data/data/hifi/mapped/CHM13/${SAMPLE_CHILD}.CHM13.haplotagged.bam \
      ${PRJ_DIR}/hifi/${SAMPLE_CHILD}.CHM13.haplotagged.bam

    aws s3 cp --no-sign-request \
      s3://platinum-pedigree-data/data/hifi/mapped/CHM13/${SAMPLE_CHILD}.CHM13.haplotagged.bam.bai \
      ${PRJ_DIR}/hifi/${SAMPLE_CHILD}.CHM13.haplotagged.bam.bai


    ## CHILD - REF GRch38
    # aws s3 cp --no-sign-request \
    #   s3://platinum-pedigree-data/data/hifi/mapped/GRCh38/${SAMPLE_CHILD}.GRCh38.haplotagged.bam \
    #   ${PRJ_DIR}/hifi/${SAMPLE_CHILD}.GRCh38.haplotagged.bam


    # aws s3 cp --no-sign-request \
    #   s3://platinum-pedigree-data/data/hifi/mapped/GRCh38/${SAMPLE_CHILD}.GRCh38.haplotagged.bam.bai \
    #   ${PRJ_DIR}/hifi/${SAMPLE_CHILD}.GRCh38.haplotagged.bam.bai

    # ## Current not used
    ## =============================================================================================
    # # ONT BAMs
    # aws s3 cp --no-sign-request \
    #   s3://platinum-pedigree-data/data/ont/mapped/CHM13/${SAMPLE_CHILD}.minimap2.bam \
    #   ./${SAMPLE_CHILD}.CHM13.minimap2.bam

    # aws s3 cp --no-sign-request \
    #   s3://platinum-pedigree-data/data/ont/mapped/CHM13/${SAMPLE_CHILD}.minimap2.bam.bai \
    #   ./${SAMPLE_CHILD}.CHM13.minimap2.bam.bai

    ## =============================================================================================
    ## VCF 
    # Pedigree consistent merged small variant calls (truthset) - Not used
    # aws s3 cp --no-sign-request  \
    #   s3://platinum-pedigree-data/variants/small_variant_truthset/GRCh38/CEPH1463.GRCh38.family-truthset.ov.vcf.gz \
    #   ./small_variant_truthset/CEPH1463.GRCh38.family-truthset.ov.vcf.gz

    # aws s3 cp --no-sign-request  \
    #   s3://platinum-pedigree-data/variants/small_variant_truthset/GRCh38/hq_regions_final.bed.gz \
    #   ./small_variant_truthset/CEPH1463.GRCh38.family-truthset.ov.vcf.gz
    
    # Dragen (Illumina) calls 
    # aws s3 cp --no-sign-request  \
    # s3://platinum-pedigree-data/variants/small_variants/CHM13/CEPH1463.CHM13.illumina-dragen.oa.vcf.gz \
    # ${ILLUM_DIR}/CEPH1463.CHM13.illumina-dragen.oa.vcf.gz

    # aws s3 cp --no-sign-request \
    # s3://platinum-pedigree-data/variants/small_variants/GRCh38/CEPH1463.GRCh38.illumina-dragen.oa.vcf.gz \
    # ${ILLUM_DIR}/CEPH1463.GRCh38.illumina-dragen.oa.vcf.gz

}


download_hifi_data_job() {

    local OUT=./download_hifi_${SAMPLE_CHILD}.slurm
    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J download_hifi_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH -A r00379

set -euo pipefail

module load aws-cli/2.25.5

SAMPLE=${SAMPLE_CHILD}
PRJ_DIR=${PRJ_DIR}

OUT_DIR=\${PRJ_DIR}/hifi
mkdir -p \${OUT_DIR}

echo "[download] Downloading HiFi BAM for \${SAMPLE}"

aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/data/hifi/mapped/CHM13/\${SAMPLE}.CHM13.haplotagged.bam \
  \${OUT_DIR}/\${SAMPLE}.CHM13.haplotagged.bam

aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/data/hifi/mapped/CHM13/\${SAMPLE}.CHM13.haplotagged.bam.bai \
  \${OUT_DIR}/\${SAMPLE}.CHM13.haplotagged.bam.bai

echo "[download] Done"
EOF

    chmod +x "${OUT}"

    ## Extract the id of the job

    # JOBID=$(sbatch --parsable "${OUT}")
    sbatch "${OUT}"

    # Communicate it with next process
    # echo "${JOBID}"

    rm "${OUT}"
}



#####################################################3
## Preprocess vcf files
#######################################################


# optional
# ------------------------------------------------------------------------------
# # Check disk memory used, one level deep, then sort by size
# du -sh * .[^.]* 2>/dev/null | sort -h

# # Check file in aws platinum pedigree
# aws s3 ls --no-sign-request s3://platinum-pedigree-data/
# ------------------------------------------------------------------------------

# # link to reference
# https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz


##############################################
# 2. Preprocess VCFs (split mom / child vcfs)
##############################################
split_vcfs() {

    local SAMPLE=$1
    local REFGENOME=$2
    local infile=$3

    bcftools view -s "${SAMPLE}" \
      -Oz -o "${ILLUM_DIR}/${SAMPLE}.${REFGENOME}.illumina-dragen.oa.vcf.gz" \
      "${infile}"

    tabix -p vcf "${ILLUM_DIR}/${SAMPLE}.${REFGENOME}.illumina-dragen.oa.vcf.gz"
}



#####################################################3
## 3. Preprocess vcf files (remove PS, phasing artifacts)
#######################################################
clean_original_vcf_illumina() {

    local SAMPLE=$1
    local REFGENOME=$2


    local INFILE=${ILLUM_DIR}/${SAMPLE}.${REFGENOME}.illumina-dragen.oa.vcf.gz
    local OUTFILE=${ILLUM_DIR}/${SAMPLE}.${REFGENOME}.illumina-dragen.oa.unphased.noPS.vcf.gz

    echo "[clean_original_vcf] Cleaning ${SAMPLE} (ref file : ${REF})"

    bcftools +setGT "${INFILE}" -Ou -- -t a -n u \
    | bcftools annotate -x FORMAT/PS \
    | bcftools view -Oz -o "${OUTFILE}"

    bcftools index "${OUTFILE}"
}


###################################################
## 3.5 (optional) Preprocess vcfs ( normalize and compare illumina and hifi vcfs)
###################################################
## These output files are for analysis of False positive. 
## They don't have to be used for phasing as Whatshap will realign things
norm_and_compare_vcfs(){

    local ILLUM_VCF_CHILD_in=${ILLUM_DIR}/${SAMPLE_CHILD}.CHM13.illumina-dragen.oa.unphased.noPS.vcf.gz
    local ILLUM_VCF_CHILD_out=${ILLUM_DIR}/${SAMPLE_CHILD}.CHM13.illumina-dragen.oa.unphased.noPS.normed.vcf.gz

    local HIFI_VCF_CHILD_in=${HIFI_DIR}/${SAMPLE_CHILD}.CHM13.deepvariant.glnexus.oa.vcf.gz
    local HIFI_VCF_CHILD_out=${HIFI_DIR}/${SAMPLE_CHILD}.CHM13.deepvariant.glnexus.oa.normed.vcf.gz

    local CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

    echo "Normalizing Illumina VCF (autosomes + X,Y)"
    bcftools norm \
      -f ${REF} \
      -m -any \
      -r ${CHRS} \
      ${ILLUM_VCF_CHILD_in} \
      -Oz -o ${ILLUM_VCF_CHILD_out}

    echo "Normalizing HiFi VCF (autosomes + X,Y)"
    bcftools norm \
      -f ${REF} \
      -m -any \
      -r ${CHRS} \
      ${HIFI_VCF_CHILD_in} \
      -Oz -o ${HIFI_VCF_CHILD_out}

    bcftools index ${ILLUM_VCF_CHILD_out}
    bcftools index ${HIFI_VCF_CHILD_out}
}


##############################################
# 4. Phase the child's illumina vcf using long reads (Slurm job)
##############################################
generate_phasing_job() {
    local OUT=./aws-data/variants/small_variants/build_hapl_${SAMPLE_CHILD}.slurm

    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J build_hapl_${SAMPLE_CHILD}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME}
#SBATCH --mem=8G
#SBATCH -A r00379

set -euo pipefail

REF=${PRJ_DIR}/reference/chm13v2.0_maskedY_rCRS.fa
SAMPLE=${SAMPLE_CHILD}
ILLUM_VCF=${ILLUM_DIR}/\${SAMPLE}.CHM13.illumina-dragen.oa.unphased.noPS.vcf.gz
HIFI_BAM=${PRJ_DIR}/hifi/\${SAMPLE}.CHM13.haplotagged.bam
module load conda
conda activate whatshap-env

whatshap phase \
  --reference \${REF} \
  -o \${SAMPLE}.illumVCF_hifiOnly.phased.vcf.gz \
  \${ILLUM_VCF} \${HIFI_BAM}

whatshap stats \
  --tsv=\${SAMPLE}.illumVCF_hifiOnly.stats.tsv \
  --block-list=\${SAMPLE}.illumVCF_hifiOnly.blocks.tsv \
  \${SAMPLE}.illumVCF_hifiOnly.phased.vcf.gz

echo "Done"
EOF
    # cat "${OUT}"
    chmod +x "${OUT}"
    
    # if [[ -n "${DEPENDENCY_JOBID}" ]]; then
    #     JOBID=$(sbatch --parsable \
    #         --dependency=afterok:${DEPENDENCY_JOBID} \
    #         "${OUT}")
    # else
    #     JOBID=$(sbatch --parsable "${OUT}")
    # fi
    sbatch "${OUT}"
    rm "${OUT}"
    echo "${JOBID}"

}

##############################################
# 5. Analyze phased mom/child variants output files
##############################################


merge_unphased-parent_phased-child_vcfs() {
    # clean up
    rm -f ./slurm*out
    rm -f ./typescript

    # Make and move phased vcfs to correct directories
    local PHASED_VCF=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged

    mkdir -p ${PHASED_VCF}
    mkdir -p ${MERGED_PHASED_VCF}

    # move child vcf to correct place
    mv ./${SAMPLE_CHILD}.illumVCF_hifiOnly* ${PHASED_VCF}/
    
    # Index samples before merging
    bcftools index -f ${ILLUM_DIR}/${SAMPLE_PARENT}.CHM13.illumina-dragen.oa.unphased.noPS.vcf.gz
    bcftools index -f ${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_hifiOnly.phased.vcf.gz

    # Merge samples into one vcf
    echo "merging samples"
    bcftools merge -m none -Oz \
        -o ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz \
        ${ILLUM_DIR}/${SAMPLE_PARENT}.CHM13.illumina-dragen.oa.unphased.noPS.vcf.gz \
        ${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_hifiOnly.phased.vcf.gz

    bcftools index ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.vcf.gz #index
    echo "done"
}



extract_phased_snp() {
    local PHASED_VCF=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf
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

    python ${PRJ_DIR}/variants/src/fix_PhaseSet.py \
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
  local PHASED_VCF=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf
  local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks
  local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis

  mkdir -p ${MISMATCH_ANALYSIS}
  
  echo " Counting alleles per haplotype block"
  conda activate whatshap-env
  python ${PRJ_DIR}/variants/src/count_mismatches.py \
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
    local PHASED_VCF=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf
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
  
  local HIFI_BAM=${PRJ_DIR}/hifi/${SAMPLE_CHILD}.CHM13.haplotagged.bam
  local DNM_BED=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.bed

  python ${PRJ_DIR}/variants/src/dnmc_readcheck.py \
  ${SAMPLE_CHILD} \
  ${HIFI_BAM} \
  ${DNM_BED} \
  ${MIN_BASE_QUAL} \
  ${MIN_MAP_QUAL} \
  ${WINDW} ${ALT_READ_COUNT} F
}



### ===================================================================================
## Second validation step
## ======================================================================================


## =====Remove later========
## check all the candidates suggested

# samtools mpileup -l ../variants/small_variants/NA12879_phasedvcf/mismatch_analysis/NA12878_NA12879_dnmc.bed \
# NA12879.CHM13.haplotagged.bam -f ../reference/chm13v2.0_maskedY_rCRS.fa 
## ======================================================================



##############################################
regenerate_phasing_job() {
    # local DEPENDENCY_JOBID=$1
    local OUT=./aws-data/variants/small_variants/build_hapl_${SAMPLE_CHILD}.slurm
    local DNM_BED=${PRJ_DIR}/variants/small_variants/NA12879_hifi_child_only_filtered_dnmc.bed

    cat << EOF > "${OUT}"
#!/bin/bash
#SBATCH -J build_hapl_${SAMPLE_CHILD}_localphase
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME}
#SBATCH --mem=8G
#SBATCH -A r00379

set -euo pipefail

module load bcftools
module load conda
conda activate whatshap-env


REF=${PRJ_DIR}/reference/chm13v2.0_maskedY_rCRS.fa
SAMPLE=${SAMPLE_CHILD}
v=${v}

ILLUM_VCF=${ILLUM_DIR}/\${SAMPLE}.CHM13.illumina-dragen.oa.unphased.noPS.vcf.gz
HIFI_BAM=${PRJ_DIR}/hifi/\${SAMPLE}.CHM13.haplotagged.bam

REGIONS_BED=\${SAMPLE}_dnm_plusminus\${v}kb.bed
LOCAL_VCF=\${SAMPLE}_dnm_plusminus\${v}kb.vcf.gz

echo "[local phase] Subsetting Illumina VCF to local regions"

bcftools view \
  -R \${REGIONS_BED} \
  -Oz -o \${LOCAL_VCF} \
  \${ILLUM_VCF}

bcftools index \${LOCAL_VCF}

echo "[local phase] Running whatshap on local regions only"

whatshap phase \
  --reference \${REF} \
  -o \${SAMPLE}.illumVCF_hifiOnly.phased.\${v}kb.vcf.gz \
  \${LOCAL_VCF} \${HIFI_BAM}

whatshap stats \
  --tsv=\${SAMPLE}.illumVCF_hifiOnly.stats.\${v}kb.tsv \
  --block-list=\${SAMPLE}.illumVCF_hifiOnly.blocks.\${v}kb.tsv \
  \${SAMPLE}.illumVCF_hifiOnly.phased.\${v}kb.vcf.gz

echo "Done local DNM re-phasing"
EOF

    chmod +x "${OUT}"
    # if [[ -n "${DEPENDENCY_JOBID}" ]]; then
    #     JOBID=$(sbatch --parsable \
    #         --dependency=afterok:${DEPENDENCY_JOBID} \
    #         "${OUT}")
    # else
    #     JOBID=$(sbatch --parsable "${OUT}")
    # fi
    sbatch "${OUT}"
    rm "${OUT}"
}




remerge_unphased-parent_phased-child_vcfs() {
    # clean up
    rm -f ./slurm*out
    rm -f ./typescript

    # Make and move phased vcfs to correct directories
    local PHASED_VCF=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged


    mkdir -p ${PHASED_VCF}
    mkdir -p ${MERGED_PHASED_VCF}

    # move child vcf to correct place
    mv -f ./${SAMPLE_CHILD}.illumVCF_hifiOnly* ${PHASED_VCF}/
    
    # Index samples before merging
    bcftools index -f ${ILLUM_DIR}/${SAMPLE_PARENT}.CHM13.illumina-dragen.oa.unphased.noPS.vcf.gz
    bcftools index -f ${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_hifiOnly.phased.${v}kb.vcf.gz

    # Merge samples into one vcf
    echo "merging samples"
    bcftools merge -m none -Oz \
        -o ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.${v}kb.vcf.gz \
        ${ILLUM_DIR}/${SAMPLE_PARENT}.CHM13.illumina-dragen.oa.unphased.noPS.vcf.gz \
        ${PHASED_VCF}/${SAMPLE_CHILD}.illumVCF_hifiOnly.phased.${v}kb.vcf.gz

    bcftools index ${MERGED_PHASED_VCF}/${SAMPLE_PARENT}_${SAMPLE_CHILD}.merged.${v}kb.vcf.gz #index
    echo "done"
}


reextract_phased_snp() {
    local PHASED_VCF=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf
    local MERGED_PHASED_VCF=${PHASED_VCF}/merged
    local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks
    

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
    python ${PRJ_DIR}/variants/src/fix_PhaseSet.py \
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
  local PHASED_VCF=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf
  local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks
  local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis
  local RECOUNT="T"
  local DNMC_FILE=${PHASED_VCF}/mismatch_analysis/${SAMPLE_PARENT}_${SAMPLE_CHILD}_dnmc.tsv
  
  conda activate whatshap-env
  python "${PRJ_DIR}/variants/src/count_mismatches.py" \
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
  local SLICE_DIR=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis/sliced/
  
  mkdir -p ${SLICE_DIR}

  mv ${SAMPLE_CHILD}*plusminus* ${SLICE_DIR}
  mv ${SAMPLE_CHILD}*LR_validated* ${SLICE_DIR}

}


calculate_callable_genome() {
  local PHASED_VCF=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf
  local EXTRACT_HAPLBLOCK=${PHASED_VCF}/HTblocks
  local MISMATCH_ANALYSIS=${PHASED_VCF}/mismatch_analysis
  local POST_ANALYSIS=${PHASED_VCF}/mismatch_analysis/denum_calcul
  
  mkdir -p ${POST_ANALYSIS}
  conda activate whatshap-env
  python "${PRJ_DIR}/variants/src/callable_genome.py" \
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
    local PHASED_VCF=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf
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
  
  local HIFI_BAM=${PRJ_DIR}/hifi/${SAMPLE_CHILD}.CHM13.haplotagged.bam
  local DNM_BED=${PRJ_DIR}/variants/small_variants/${SAMPLE_CHILD}_phasedvcf/mismatch_analysis/denum_calcul/${SAMPLE_PARENT}_${SAMPLE_CHILD}_hetc.bed

  python ${PRJ_DIR}/variants/src/dnmc_readcheck.py \
  ${SAMPLE_CHILD} \
  ${HIFI_BAM} \
  ${DNM_BED} \
  ${MIN_BASE_QUAL} \
  ${MIN_MAP_QUAL} \
  ${WINDW} ${ALT_READ_COUNT} ${VERBOSE}
}


##################################################
# Main
##################################################


main() {

    load_modules  # always

    for PART in "${PARTS[@]}"; do

        ############################################
        # PART 1 — Data Preparation
        ############################################
        if [[ "$PART" == "1" ]]; then
            echo "========== PART 1: Data preparation =========="

            # Download data from AWS
            download_hifi_data_job

            # Create a vcf for each individual
            split_vcfs ${SAMPLE_PARENT} CHM13 ${ILLUM_DIR}/CEPH1463.CHM13.illumina-dragen.oa.vcf.gz
            split_vcfs ${SAMPLE_CHILD} CHM13 ${ILLUM_DIR}/CEPH1463.CHM13.illumina-dragen.oa.vcf.gz

            # Remove "PS" tags or "|" in genotype
            clean_original_vcf_illumina ${SAMPLE_CHILD} CHM13
            clean_original_vcf_illumina ${SAMPLE_PARENT} CHM13

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
            merge_unphased-parent_phased-child_vcfs

            # Extract only positions where there's a Phase Set (PS) associated
            extract_phased_snp

            # Accounting on PS blocks
            count_shared_alleles_per_PS_block

            # Create a list of dnm candidates
            filter_dnm_candidates

            # Additional filters with hifi suppport
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
        fi


    done

    echo "Pipeline finished."
}