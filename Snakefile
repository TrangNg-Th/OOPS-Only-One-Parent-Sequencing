"""
OOPS – Only One Parent Sequencing
Snakemake pipeline for BAM/VCF analysis using samtools, bcftools, and whatshap.

Steps
-----
1. sort    – Sort BAM by coordinate (samtools sort)
2. index   – Index sorted BAM (samtools index)
3. stats   – Alignment statistics (samtools flagstat)
4. call    – Variant calling (bcftools mpileup | bcftools call)
5. filter  – Hard-filter variants by QUAL and depth (bcftools filter/view)
6. phase   – Phase variants per sample (whatshap phase)
7. haplotag– Tag reads in BAM with haplotype assignments (whatshap haplotag)
8. phase_stats – Summary statistics for phased VCF (whatshap stats)
"""

configfile: "config/config.yml"

SAMPLES   = config["samples"]
REF       = config["reference"]
BAM_DIR   = config["bam_dir"]
RESULTS   = config["results_dir"]

# ---------------------------------------------------------------------------
# Top-level target
# ---------------------------------------------------------------------------
rule all:
    input:
        expand("{results}/flagstat/{sample}.flagstat.txt",  results=RESULTS, sample=SAMPLES),
        expand("{results}/vcf/{sample}.filtered.vcf.gz",    results=RESULTS, sample=SAMPLES),
        expand("{results}/phased/{sample}.phased.vcf.gz",   results=RESULTS, sample=SAMPLES),
        expand("{results}/haplotagged/{sample}.haplotagged.bam", results=RESULTS, sample=SAMPLES),
        expand("{results}/phase_stats/{sample}.stats.tsv",  results=RESULTS, sample=SAMPLES),


# ---------------------------------------------------------------------------
# 1. Sort BAM
# ---------------------------------------------------------------------------
rule sort_bam:
    input:
        bam = BAM_DIR + "/{sample}.bam",
    output:
        bam = RESULTS + "/sorted/{sample}.sorted.bam",
    log:
        RESULTS + "/logs/sort/{sample}.log",
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} 2>{log}"


# ---------------------------------------------------------------------------
# 2. Index sorted BAM
# ---------------------------------------------------------------------------
rule index_bam:
    input:
        bam = RESULTS + "/sorted/{sample}.sorted.bam",
    output:
        bai = RESULTS + "/sorted/{sample}.sorted.bam.bai",
    log:
        RESULTS + "/logs/index/{sample}.log",
    shell:
        "samtools index {input.bam} 2>{log}"


# ---------------------------------------------------------------------------
# 3. Alignment statistics
# ---------------------------------------------------------------------------
rule flagstat:
    input:
        bam = RESULTS + "/sorted/{sample}.sorted.bam",
        bai = RESULTS + "/sorted/{sample}.sorted.bam.bai",
    output:
        txt = RESULTS + "/flagstat/{sample}.flagstat.txt",
    log:
        RESULTS + "/logs/flagstat/{sample}.log",
    shell:
        "samtools flagstat {input.bam} > {output.txt} 2>{log}"


# ---------------------------------------------------------------------------
# 4. Variant calling
# ---------------------------------------------------------------------------
rule call_variants:
    input:
        bam = RESULTS + "/sorted/{sample}.sorted.bam",
        bai = RESULTS + "/sorted/{sample}.sorted.bam.bai",
        ref = REF,
    output:
        vcf = RESULTS + "/vcf/{sample}.raw.vcf.gz",
    log:
        RESULTS + "/logs/call/{sample}.log",
    params:
        min_bq = config["min_base_quality"],
        min_mq = config["min_map_quality"],
        model  = config["calling_model"],
    threads: 4
    shell:
        """
        bcftools mpileup \
            --fasta-ref {input.ref} \
            --min-BQ {params.min_bq} \
            --min-MQ {params.min_mq} \
            --output-type u \
            {input.bam} \
        | bcftools call \
            -{params.model} \
            --variants-only \
            --output-type z \
            --output {output.vcf} \
            2>{log}
        bcftools index --tbi {output.vcf}
        """


# ---------------------------------------------------------------------------
# 5. Filter variants
# ---------------------------------------------------------------------------
rule filter_variants:
    input:
        vcf = RESULTS + "/vcf/{sample}.raw.vcf.gz",
    output:
        vcf = RESULTS + "/vcf/{sample}.filtered.vcf.gz",
    log:
        RESULTS + "/logs/filter/{sample}.log",
    params:
        min_qual  = config["min_variant_quality"],
        min_depth = config["min_depth"],
    shell:
        """
        bcftools filter \
            --exclude 'QUAL<{params.min_qual} || INFO/DP<{params.min_depth}' \
            --output-type z \
            --output {output.vcf} \
            {input.vcf} \
            2>{log}
        bcftools index --tbi {output.vcf}
        """


# ---------------------------------------------------------------------------
# 6. Phase variants with whatshap
# ---------------------------------------------------------------------------
rule phase_variants:
    input:
        vcf = RESULTS + "/vcf/{sample}.filtered.vcf.gz",
        bam = RESULTS + "/sorted/{sample}.sorted.bam",
        bai = RESULTS + "/sorted/{sample}.sorted.bam.bai",
        ref = REF,
    output:
        vcf = RESULTS + "/phased/{sample}.phased.vcf.gz",
    log:
        RESULTS + "/logs/phase/{sample}.log",
    params:
        ploidy   = config["ploidy"],
        max_cov  = config["whatshap_max_coverage"],
    shell:
        """
        whatshap phase \
            --reference {input.ref} \
            --ploidy {params.ploidy} \
            --max-coverage {params.max_cov} \
            --output {output.vcf} \
            {input.vcf} \
            {input.bam} \
            2>{log}
        bcftools index --tbi {output.vcf}
        """


# ---------------------------------------------------------------------------
# 7. Haplotag reads
# ---------------------------------------------------------------------------
rule haplotag:
    input:
        vcf = RESULTS + "/phased/{sample}.phased.vcf.gz",
        bam = RESULTS + "/sorted/{sample}.sorted.bam",
        bai = RESULTS + "/sorted/{sample}.sorted.bam.bai",
        ref = REF,
    output:
        bam = RESULTS + "/haplotagged/{sample}.haplotagged.bam",
    log:
        RESULTS + "/logs/haplotag/{sample}.log",
    shell:
        """
        whatshap haplotag \
            --reference {input.ref} \
            --output {output.bam} \
            {input.vcf} \
            {input.bam} \
            2>{log}
        samtools index {output.bam}
        """


# ---------------------------------------------------------------------------
# 8. Phasing statistics
# ---------------------------------------------------------------------------
rule phase_stats:
    input:
        vcf = RESULTS + "/phased/{sample}.phased.vcf.gz",
    output:
        tsv = RESULTS + "/phase_stats/{sample}.stats.tsv",
    log:
        RESULTS + "/logs/phase_stats/{sample}.log",
    shell:
        "whatshap stats --tsv {output.tsv} {input.vcf} 2>{log}"
