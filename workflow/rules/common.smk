import pandas as pd
import os
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"   # config オブジェクトを作成。config.yamlが辞書形式で格納される。
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)
validate(units, schema="../schemas/units.schema.yaml")

##### Wildcard constraints #####
# wildcard_constraints: でglobalに定義。
# '|' は正規表現をORで繋ぐ表現方法
# | でつないだ単語のみが扱う範囲となる。例えば、vartypeは正規表現の(snvs | indels)と等価。
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index)


##### Helper functions #####

def get_fasta_basename():
    return os.path.basename(config["ref"]["genome"])


def get_fastq(wildcards):
    """
    Get fastq files of given sample-unit.
    # fastq fileのパスをreturen 
    """
    fastqs = units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()   # dropna()でNaNを含む行を削除
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def get_fastq_for_fastp(wildcards):
    """
    Get fastq files of given sample-unit.
    # fastq fileのパスをreturen 
    """
    fastqs = units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()   # dropna()でNaNを含む行を削除
    if len(fastqs) == 2:
        return {"sample": [fastqs.fq1, fastqs.fq2]}


def is_single_end(sample):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("data/output/trimmed/{sample}.{group}.fq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "data/output/trimmed/{sample}.fq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
                #   "recal/{sample}.bam",
                  "data/output/dedup/{sample}.bam",
                  sample=wildcards.sample)


def get_sample_bams_bamindex(wildcards):
    # rule call_variantsにおいてunpack
    # baiもrule call_variantsの入力に入れておかなければrule samtools indexが呼び出されない。
    """Get all aligned reads of given sample."""
    return {"bam": expand(
                  "data/output/dedup/{sample}.bam",
                  sample=wildcards.sample),
            "bai": expand(
                  "data/output/dedup/{sample}.bam.bai",
                  sample=wildcards.sample)
                  }
