if config["processing"]["trimming"]=="trimmomatic":
    rule trim_reads_se:
        input:
            unpack(get_fastq)
        output:
            temp("data/output/trimmed/{sample}.fq.gz")
        params:
            extra="",
            **config["params"]["trimmomatic"]["se"]
        log:
            "logs/trimmomatic/{sample}.log"
        wrapper:
            "0.50.4/bio/trimmomatic/se"


    rule trim_reads_pe:
        input:
            unpack(get_fastq)
        output:
            r1=temp("data/output/trimmed/{sample}.1.fq.gz"),
            r2=temp("data/output/trimmed/{sample}.2.fq.gz"),
            r1_unpaired=temp("data/output/trimmed/{sample}.1.unpaired.fq.gz"),
            r2_unpaired=temp("data/output/trimmed/{sample}.2.unpaired.fq.gz"),
            trimlog="data/output/trimmed/{sample}.trimlog.txt"
        params:
            extra=lambda w, output: "-trimlog {}".format(output.trimlog),
            **config["params"]["trimmomatic"]["pe"]
        log:
            "logs/trimmomatic/{sample}.log"
        wrapper:
            "0.50.4/bio/trimmomatic/pe"


if config["processing"]["trimming"]=="fastp":
    rule fastp_pe: 
        input: 
            unpack(get_fastq_for_fastp)
        output:
            trimmed=temp(["data/output/trimmed/{sample}.1.fq.gz", "data/output/trimmed/{sample}.2.fq.gz"]),
            html="data/output/qc/fastp/{sample}_fastp.html",
            json="data/output/qc/fastp/{sample}_fastp.json"
        log:
            "logs/fastp/{sample}.log"
        params:
            extra=" ".join([*config["params"]["fastp"]])
        threads: 2
        wrapper:
            "0.50.4/bio/fastp"


rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output
    output:
        temp("data/output/mapped/{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params: # lambda 参考：https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-3-input-functions
        index=config["ref"]["genome"],
        extra=lambda wildcards: get_read_group(wildcards)\
         + " " + " ".join([*config["params"]["bwa"]["mem"]]),
        sort="samtools",
        sort_order="coordinate"
    threads: 2
    wrapper:
        "0.50.4/bio/bwa/mem"


rule mark_duplicates:
    input:
        "data/output/mapped/{sample}.sorted.bam"
    output:
        bam=protected("data/output/dedup/{sample}.bam"),
        metrics="data/output/qc/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.50.4/bio/picard/markduplicates"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        protected("{prefix}.bam.bai")
    wrapper:
        "0.50.4/bio/samtools/index"


rule RealignerTargetCreator:
    input:
        unpack(get_sample_bams_bamindex),
        ref = config["ref"]["genome"],
        idx=rules.genome_dict.output
    output:
        target_list = protected("data/output/realign/{sample}.list")
    log:
        "logs/gatk/realign/{sample}_target.log"
    benchmark:
        "data/output/benchmark/gatk/realing/{sample}_target.tsv"
    threads: 2
    conda:
        "../envs/gatklite.yaml"
    shell:
        "java -jar $CONDA_DEFAULT_ENV/jar/GenomeAnalysisTKLite.jar "
        "-T RealignerTargetCreator "
        "-nt {threads} "
        "--validation_strictness LENIENT "
        "-R {input.ref} "
        "-I {input.bam} "
        "-o {output.target_list} | tee -a {log}"


rule IndelRealigner:
    input:
        unpack(get_sample_bams_bamindex),
        target_list = "data/output/realign/{sample}.list",
        ref = config["ref"]["genome"],
    output:
        bam = protected("data/output/realign/{sample}_align.bam"),
        bai = protected("data/output/realign/{sample}_align.bai"),
    log:
        "logs/gatk/realign/{sample}.log"
    benchmark:
        "data/output/benchmark/gatk/realing/{sample}_realign.tsv"
    conda:
        "../envs/gatklite.yaml"
    shell:
        "java -jar $CONDA_DEFAULT_ENV/jar/GenomeAnalysisTKLite.jar "
        "-T IndelRealigner "
        "--validation_strictness LENIENT "
        "--read_filter NotPrimaryAlignment "
        "-R {input.ref} "
        "-I {input.bam} "
        "-targetIntervals {input.target_list} "
        "-o {output.bam} | tee -a {log}"