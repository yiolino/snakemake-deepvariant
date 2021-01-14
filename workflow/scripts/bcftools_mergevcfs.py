from snakemake.shell import shell

inputs = " ".join(f for f in snakemake.input.vcfs)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

shell(
    "bcftools"
    " merge"
    " {inputs}"
    " | pigz -c"
    " > {snakemake.output[0]}"
    " {log}"
)