__author__ = "Tetsuro Hisayoshi"
__copyright__ = "Copyright 2020, Tetsuro Hisayoshi"
__email__ = "hisayoshi0530@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

log_dir = os.path.dirname(str(snakemake.log))
output_dir = os.path.dirname(str(snakemake.output[0]))

# sample basename
basename = os.path.splitext(os.path.basename(str(snakemake.input.bam)))[0]
print(basename)

split_inputs = " ".join(str(x) for x in range(0, int(snakemake.threads)))


shell(
    "(BIN_DIR=$(ls -d $CONDA_DEFAULT_ENV/share/deepvariant*/binaries/Deepvariant/*/DeepVariant*) \n"
    "mkdir -p {output_dir}/{basename} \n"
    "parallel --eta --halt 2 --joblog {log_dir}/log --res {log_dir} "
    "python $BIN_DIR/make_examples.zip "
    "--mode calling --ref {snakemake.input.ref} --reads {snakemake.input.bam} "
    "--examples {output_dir}/{basename}/{basename}.tfrecord@{snakemake.threads}.gz "
    "--gvcf {output_dir}/{basename}/{basename}.gvcf.tfrecord@{snakemake.threads}.gz "
    "--task {{}} "
    "::: {split_inputs} \n"
    "dv_call_variants.py "
    "--cores {snakemake.threads} "
    "--outfile {output_dir}/{basename}/{basename}.tmp "
    "--sample {basename} "
    "--examples {output_dir}/{basename} "
    "--model {snakemake.params.model} \n"
    "python $BIN_DIR/postprocess_variants.zip "
    "--ref {snakemake.input.ref} "
    "--infile {output_dir}/{basename}/{basename}.tmp "
    "--outfile {snakemake.output.vcf} "
    "--nonvariant_site_tfrecord_path {output_dir}/{basename}/{basename}.gvcf.tfrecord@{snakemake.threads}.gz "
    "--gvcf_outfile {snakemake.output.gvcf} \n"
    "rm -rf {output_dir}/{basename} ) {log}"
)