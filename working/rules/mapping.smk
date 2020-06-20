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
