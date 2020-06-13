# snakemake-deepvariant
This is sankemake-wrapper of variant calling workflow using [DeepVariant](https://github.com/google/deepvariant).  

# Authors
tetsuro90

# How to use

#### Step 1: Obtain a copy of this workflow
* Create a new github repository using this workflow as a template.
* Clone the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 2: Move your fastq file to `data/reads/`

#### Step 3: Move your fasta file to `data/ref/`

#### Step 4: edit `working/unit.tsv` and `working/sample.tsv`
An example is shown below.  
```
sample	platform	fq1	fq2
SRR1770413	ILLUMINA	data/reads/SRR1770413_1.fastq.gz	data/reads/SRR1770413_2.fastq.gz
SRR341549	ILLUMINA	data/reads/SRR341549_1.fastq.gz	data/reads/SRR341549_2.fastq.gz
```

```
sample
SRR1770413
SRR341549
```

#### Step5: create fasta index and dict
```
$ cd data/ref/
$ samtools faidx {your ref}.fa
$ picard CreateSequenceDictionary R={your ref}.fa O={your ref}.dict
```

#### Step6: run sankemake workflow
```
cd working/
snakemake -j {#CORE} --use-conda all
```

#### Step7: make a report html
Once the workflow has been successfully completed, You can output a report file.
```
snakemake --report
```

# License
MIT License

# Acknowledgements
* [DeepVariant](https://github.com/google/deepvariant). Copyright 2017 Google LLC. BSD-3-Clause license.
* [DeepVariant bioconda package](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/deepvariant). Copyright 2018 Brad Chapman. MIT license.
* [Snakemake workflow: dna-seq-gatk-variant-calling](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling). Copyright 2018, Johannes Köster. MIT license.
* [fastp](https://github.com/OpenGene/fastp). Copyright 2017 OpenGene - Open Source Genetics Toolbox. MIT license.
* [samtools](https://github.com/samtools/samtools). Copyright 2008-2019 Genome Research Ltd. The MIT/Expat License.
* [gatk4](https://github.com/broadinstitute/gatk). Copyright 2009-2020, Broad Institute, Inc. BSD-3-Clause license.
* [bwa](https://github.com/lh3/bwa). Copyright 2009, Li H. and Durbin R. GPLv3 License.
 
