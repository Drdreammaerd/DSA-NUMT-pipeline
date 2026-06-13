# ============================================================
# NUMT DSA Pipeline Environment
# ============================================================
# Stage 1: Extract LAST tools from last-train
FROM --platform=linux/amd64 dreammaerd/last-train:v5 AS last-builder

# Stage 2: Final clean image using micromamba to prevent OOM
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.6

USER root
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago
ENV MAMBA_ROOT_PREFIX=/opt/conda

# 1. Install prerequisites
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# 2. Copy LAST and custom scripts from last-builder
COPY --from=last-builder /usr/local/bin/last* /usr/local/bin/
COPY --from=last-builder /usr/local/bin/maf* /usr/local/bin/
COPY --from=last-builder /usr/local/bin/seg* /usr/local/bin/
COPY --from=last-builder /usr/local/bin/merge /usr/local/bin/
COPY --from=last-builder /usr/local/bin/remvrrna /usr/local/bin/
COPY --from=last-builder /usr/local/bin/fix-protein-strand /usr/local/bin/

# Fix syntax error in fix-protein-strand
RUN sed -i 's/if \[$# -ne argNum\]; then/if \[ $# -ne $argNum \]; then/g' /usr/local/bin/fix-protein-strand

# 3. Install Dependencies using micromamba
USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    snakemake \
    pandas \
    bedtools \
    samtools \
    pysam \
    blast \
    r-base \
    r-optparse \
    r-data.table \
    r-stringr \
    bioconductor-biostrings \
    bioconductor-rtracklayer \
    bioconductor-genomicranges \
    && micromamba clean -a -y

# 4. Install Snakemake cluster executor plugin
USER root
RUN /opt/conda/bin/pip install --no-cache-dir \
    snakemake-executor-plugin-cluster-generic

# 5. Convenient symlinks
RUN ln -sf /opt/conda/bin/python  /usr/local/bin/python  && \
    ln -sf /opt/conda/bin/python  /usr/local/bin/python3 && \
    ln -sf /opt/conda/bin/snakemake /usr/local/bin/snakemake

# 6. Set Working Directory and copy pipeline files
USER root
WORKDIR /opt/numt-dsa-pipeline

COPY scripts/   /opt/numt-dsa-pipeline/scripts/
COPY resources/ /opt/numt-dsa-pipeline/resources/
COPY helpers/   /opt/numt-dsa-pipeline/helpers/
COPY Snakefile  /opt/numt-dsa-pipeline/Snakefile

RUN chmod -R a+rX /opt/numt-dsa-pipeline/ && \
    chmod +x /opt/numt-dsa-pipeline/scripts/*.sh /opt/numt-dsa-pipeline/scripts/*.py /usr/local/bin/* 2>/dev/null || true

RUN /opt/conda/bin/makeblastdb -in /opt/numt-dsa-pipeline/resources/chrM.fa -dbtype nucl && \
    /opt/conda/bin/makeblastdb -in /opt/numt-dsa-pipeline/resources/mito_proteins.fa -dbtype prot && \
    chmod a+r /opt/numt-dsa-pipeline/resources/chrM.fa.* /opt/numt-dsa-pipeline/resources/mito_proteins.fa.*

# 7. Environment & Entrypoint
ENV PATH="/opt/conda/bin:${PATH}"

RUN printf '#!/bin/bash\neval "$(micromamba shell hook --shell bash)"\nmicromamba activate base\nexec "$@"\n' \
    > /usr/local/bin/entrypoint.sh && \
    chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["/bin/bash"]
