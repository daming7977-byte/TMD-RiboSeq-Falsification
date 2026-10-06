FROM continuumio/miniconda2:4.7.12

RUN conda install -y \
    --override-channels \
    -c bioconda \
    -c bioconda/label/cf201901 \
    -c conda-forge/label/cf201901 \
    -c defaults \
    python=2.7.15 \
    tophat=2.1.1=py27_3 \
    bowtie=1.1.2=py27_2 \
    bowtie2=2.2.5 \
 && conda clean -afy

RUN echo "============================================================" \
 && echo "FROZEN SOFTWARE CHECK" \
 && echo "============================================================" \
 && uname -m \
 && ldd --version | head -n1 \
 && python --version \
 && tophat --version \
 && bowtie --version \
 && bowtie-build --version \
 && bowtie2 --version

CMD ["/bin/bash"]
