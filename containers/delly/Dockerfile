# Native linux/arm64 (and amd64) Delly image built from pinned bioconda packages.
# Why: the osx-arm64 conda Delly 2.6.0 build dies with SIGBUS when decoding CRAM (BAM input works),
# and the published Delly images (dellytools/delly, biocontainers) are amd64-only.
FROM mambaorg/micromamba:2.9.0@sha256:e0a99b0f17a759e14c2f967dc0ca2d3a3c1ca3c62955f4d20bba770eaaf0184d
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH
