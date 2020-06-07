FROM ubuntu:18.04

ENV SNAKEMAKE_VER 5.18.0
ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends \
        # unzip \
        libfontconfig1 \
        rsync \
        tree \
        less \
        wget \
        git \
        vim \
        ssh \
        apt-transport-https \
        ca-certificates \
    && apt-get update \
    && apt-get install -y --no-install-recommends locales \
    && locale-gen en_US.UTF-8 \
    && ln -sf /usr/share/zoneinfo/Asia/Tokyo /etc/localtime \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir ${HOME}/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && conda install -c conda-forge mamba \
    && mamba create -c conda-forge -c bioconda -n snakemake snakemake==${SNAKEMAKE_VER} -y \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && echo "conda activate snakemake" >> ~/.bashrc