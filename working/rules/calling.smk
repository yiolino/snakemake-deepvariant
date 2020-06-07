rule call_variants:
    input:
        unpack(get_sample_bams_bamindex),
        ref=config["ref"]["genome"]
    output:
        vcf="data/output/deepvariant/{sample}.vcf"
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