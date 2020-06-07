rule vcf_to_tsv:
    input:
        # "data/output/annotated/all.vcf.gz"
        "data/output/filtered/all.vcf.gz"
    output:
        report("data/output/tables/calls.tsv.gz", caption="../report/calls.rst", category="Calls")
    conda:
        "../envs/rbt.yaml"
    shell:
        "bcftools view --apply-filters PASS --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}"


rule plot_stats:
    input:
        "data/output/tables/calls.tsv.gz"
    output:
        depths=report("data/output/plots/depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report("data/output/plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot-depths.py"
