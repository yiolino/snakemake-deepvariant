rule call_variants:
    input:
        unpack(get_sample_bams_bamindex),
        ref=config["ref"]["genome"]
    output:
        vcf="data/output/deepvariant/{sample}.vcf.gz",
        gvcf="data/output/deepvariant/{sample}.gvcf.gz",
    params:
        model=config["params"]["deepvariant"]["model"],
        extra=config["params"]["deepvariant"]["extra"]
    threads: 2
    log:
        "logs/deepvariant/{sample}/stdout.log"
    conda:
        "../envs/deepvariant.yaml"
    script:
        "../scripts/deepvariant.py"


rule merge_vcfs:
    input:
        vcfs=lambda w: expand("data/output/deepvariant/{sample}.vcf.gz", sample=samples.index),
    output:
        report("data/output/called/merged.vcf.gz", caption="../report/merged_vcfs.rst", category="Called SVs")
    log:
        "logs/bcftools/mergevcfs.log"
    params:
        extra=""
    conda:
        "../envs/bcftools.yaml"
    script:
        "../scripts/bcftools_mergevcfs.py"
