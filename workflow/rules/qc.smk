rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="data/output/qc/fastqc/{sample}_fastqc.html",
        zip="data/output/qc/fastqc/{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.27.1/bio/fastqc"


rule samtools_stats:
    input:
        # "recal/{sample}.bam"
        "data/output/dedup/{sample}.bam"
    output:
        "data/output/qc/samtools-stats/{sample}.txt"
    log:
        "logs/samtools-stats/{sample}.log"
    wrapper:
        "0.50.4/bio/samtools/stats"


rule multiqc:
    input:
        expand([
            "data/output/qc/fastqc/{sample}_fastqc.zip",
            "data/output/qc/samtools-stats/{sample}.txt",
            "data/output/qc/fastp/{sample}_fastp.json",
            "data/output/qc/dedup/{sample}.metrics.txt"],
            sample=list(samples.index))
    output:
        report("data/output/qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        "logs/multiqc.log"
    wrapper:
        "0.50.4/bio/multiqc"
