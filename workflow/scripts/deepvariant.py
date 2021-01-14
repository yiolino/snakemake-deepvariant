__author__ = "Tetsuro Hisayoshi"
__copyright__ = "Copyright 2020, Tetsuro Hisayoshi"
__email__ = "hisayoshi0530@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell
import tempfile

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

log_dir = os.path.dirname(snakemake.log[0])
output_dir = os.path.dirname(snakemake.output[0])

# sample basename
basename = os.path.splitext(os.path.basename(snakemake.input.bam[0]))[0]

split_inputs = " ".join(str(x) for x in range(0, int(snakemake.threads)))


with tempfile.TemporaryDirectory() as tmp_dir:
    shell(
        "(BIN_DIR=$(ls -d $CONDA_DEFAULT_ENV/share/deepvariant*/binaries/Deepvariant/*/DeepVariant*) \n"
        "parallel --eta --halt 2 --joblog {log_dir}/log --res {log_dir} "
        "python $BIN_DIR/make_examples.zip "
        "--mode calling --ref {snakemake.input.ref} --reads {snakemake.input.bam} "
        "--examples {tmp_dir}/{basename}.tfrecord@{snakemake.threads}.gz "
        "--gvcf {tmp_dir}/{basename}.gvcf.tfrecord@{snakemake.threads}.gz "
        "--task {{}} "
        "::: {split_inputs} \n"
        "dv_call_variants.py "
        "--cores {snakemake.threads} "
        "--outfile {tmp_dir}/{basename}.tmp "
        "--sample {basename} "
        "--examples {tmp_dir} "
        "--model {snakemake.params.model} \n"
        "python $BIN_DIR/postprocess_variants.zip "
        "--ref {snakemake.input.ref} "
        "--infile {tmp_dir}/{basename}.tmp "
        "--outfile {snakemake.output.vcf} "
        "--nonvariant_site_tfrecord_path {tmp_dir}/{basename}.gvcf.tfrecord@{snakemake.threads}.gz "
        "--gvcf_outfile {snakemake.output.gvcf} ) {log}"
    )
