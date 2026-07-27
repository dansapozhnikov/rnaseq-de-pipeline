# =============================================================================
# rnaseq-de-pipeline -- container image
# -----------------------------------------------------------------------------
# Base: the official Bioconductor 3.18 image (rocker-derived), which pins R 4.3,
# a Bioconductor 3.18 snapshot, a full system toolchain, and pandoc. That base
# tag is itself a reproducibility anchor: every `BiocManager::install()` below
# resolves to the 3.18 release, so the image is deterministic for a given base
# digest. (renv.lock remains the authoritative lock for host installs and is
# verified separately in CI via `renv::restore()`.)
#
# Build:  docker build -t rnaseq-de-pipeline .
# Run:    docker run --rm rnaseq-de-pipeline --validate --config config/config.example.yaml
#         docker run --rm -v "$PWD/data:/data" -v "$PWD/out:/out" rnaseq-de-pipeline \
#                    --counts /data/counts.csv --metadata /data/samples.csv \
#                    --config config/config.example.yaml --outdir /out
# =============================================================================
FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL org.opencontainers.image.title="rnaseq-de-pipeline" \
      org.opencontainers.image.description="Reproducible bulk RNA-seq differential-expression pipeline (R/DESeq2)" \
      org.opencontainers.image.source="https://github.com/dansapozhnikov/rnaseq-de-pipeline" \
      org.opencontainers.image.licenses="MIT"

WORKDIR /pipeline

# Install the analysis stack FIRST (its own layer) so editing pipeline code does
# not invalidate the expensive package-install cache. BiocManager ships in the
# base image and resolves everything against Bioconductor 3.18.
RUN R -q -e 'BiocManager::install(c( \
      "DESeq2", "apeglm", "ashr", "limma", "tximport", "SummarizedExperiment", \
      "airway", "clusterProfiler", "enrichplot", "org.Hs.eg.db", "org.Mm.eg.db", \
      "EnhancedVolcano", "matrixStats", "data.table", "readr", "ggplot2", \
      "pheatmap", "optparse", "yaml", "cli", "crayon", "logger", "rmarkdown", \
      "knitr", "testthat", "renv", "plotly", "DT", "htmlwidgets"), \
      update = FALSE, ask = FALSE, Ncpus = 4)' \
    && R -q -e 'stopifnot(requireNamespace("DESeq2"), requireNamespace("clusterProfiler"), requireNamespace("plotly"))'

# Now copy the pipeline source (cheap layer; changes here reuse the stack above).
COPY . /pipeline

# `Rscript run_pipeline.R` is the entrypoint; args are appended by `docker run`.
ENTRYPOINT ["Rscript", "run_pipeline.R"]
CMD ["--validate", "--config", "config/config.example.yaml", "--outdir", "/pipeline/results/airway_validation"]
