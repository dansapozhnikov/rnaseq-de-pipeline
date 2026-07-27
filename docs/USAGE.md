# USAGE — Reproducible Bulk RNA-seq DE Pipeline

A complete guide to installing, configuring, running, and interpreting the
pipeline. If you just want to see it work, jump to
[Quickstart](#3-running-the-pipeline) and run the `--validate` self-test.

**Contents**
1. [Overview](#1-overview)
2. [Installation](#2-installation)
3. [Running the pipeline](#3-running-the-pipeline)
4. [Inputs: count formats](#4-inputs-count-formats)
5. [The sample sheet](#5-the-sample-sheet)
6. [Configuration reference](#6-configuration-reference)
7. [CLI reference](#7-cli-reference)
8. [Outputs](#8-outputs)
9. [Reading the QC report](#9-reading-the-qc-report)
10. [Reproducibility model](#10-reproducibility-model)
11. [Troubleshooting](#11-troubleshooting)
12. [Interpreting results responsibly](#12-interpreting-results-responsibly)

---

## 1. Overview

This pipeline takes a **raw gene-level count matrix** plus a **sample sheet** and
runs the standard bulk RNA-seq DE workflow with DESeq2:

> Input → Load & align → **QC gate** → Filter → Normalize (DESeq2/VST) →
> PCA + variance-vs-metadata → Outlier & batch handling → Differential
> expression (shrinkage, thresholds) → Functional enrichment → Report

The full stage graph, including the QC decision node (a FAIL stops the run), is
in **[docs/pipeline.svg](pipeline.svg)**.

Design priorities: **correctness, reproducibility, transparency, robustness to
messy inputs, and loud QC pass/fail.**

---

## 2. Installation

Requires **R ≥ 4.2** with **Bioconductor 3.18** and **pandoc** (for the HTML
report; bundled with RStudio, or `brew install pandoc` / `apt install pandoc`).

### Option A — renv (authoritative)

`renv.lock` pins the exact package versions the pipeline was built and validated
against. From the project root:

```bash
Rscript -e 'install.packages("renv"); renv::restore()'
```

This installs everything into a project-local library (`renv/library/`), isolated
from your global R packages. First-time restore compiles/downloads the full
Bioconductor stack and can take **20–60 minutes**; subsequent restores are cached.

### Option B — conda

```bash
conda env create -f environment.yml
conda activate rnaseq-de
```

`environment.yml` mirrors the same Bioconductor 3.18 stack. `renv.lock` remains
the source of truth for exact reproducibility.

### Option C — Docker / Apptainer (most portable; recommended for HPC)

Build (or pull the CI-published image), then run the pipeline inside it — no local
R/Bioconductor install needed:

```bash
# Docker
docker build -t rnaseq-de-pipeline .
docker run --rm -v "$PWD/data:/data" -v "$PWD/out:/out" rnaseq-de-pipeline \
  --counts /data/counts.csv --metadata /data/samples.csv \
  --config config/config.example.yaml --outdir /out

# Apptainer / Singularity on HPC (no X11 needed — plots auto-use Cairo)
apptainer pull docker://ghcr.io/dansapozhnikov/rnaseq-de-pipeline:latest
apptainer run --bind "$PWD/data:/data,$PWD/out:/out" rnaseq-de-pipeline_latest.sif \
  --counts /data/counts.csv --metadata /data/samples.csv \
  --config /pipeline/config/config.example.yaml --outdir /out
```

The image is built on `bioconductor/bioconductor_docker:RELEASE_3_18`, so the
whole stack (and pandoc) is baked in and deterministic for a given base digest.

### Core packages

`DESeq2`, `apeglm` (LFC shrinkage), `limma` (batch-effect viz), `tximport`
(Salmon/kallisto), `clusterProfiler` + `org.Hs.eg.db` / `org.Mm.eg.db`
(enrichment), `ggplot2`, `pheatmap`, `EnhancedVolcano`, `plotly` + `DT`
(interactive report), `optparse`, `yaml`, `cli`/`crayon` (colored console),
`rmarkdown`/`knitr` (report).

---

## 3. Running the pipeline

### Quickstart — the airway self-test

The fastest way to confirm your install works end-to-end:

```bash
Rscript run_pipeline.R --validate --config config/config.example.yaml
```

This runs the whole pipeline on the Bioconductor `airway` dataset and asserts it
recovers the known glucocorticoid response (CRISPLD2, DUSP1, KLF15 significant
and up in treated). It exits **0** on success, **1** on any assertion failure.

### On your own data

```bash
cp config/config.example.yaml config/config.yaml   # then edit config.yaml
Rscript run_pipeline.R \
  --counts path/to/counts \
  --metadata sample_sheet.csv \
  --config config/config.yaml \
  --input-format auto \
  --outdir results/
```

A QC **FAIL** stops the run with exit code 1 (after writing the QC report). Use
`--force` to downgrade FAIL→WARN and continue (this is logged and stamped on the
report).

---

## 4. Inputs: count formats

Point `--counts` at a file or a directory. With `--input-format auto` (default)
the format is detected by inspecting the input; the chosen format and the reason
are logged. You can always force a format with `--input-format <name>`.

> **Raw counts only.** DESeq2 models raw-count overdispersion. The loaders
> **reject** non-integer/normalized input (TPM/CPM/FPKM) with a clear error —
> supply un-normalized counts (Salmon/kallisto estimated counts are the one
> exception; they are rounded to integers).

### 4.1 Plain matrix (`matrix`)

CSV/TSV with the gene ID in the first column and one numeric column per sample.

```
gene,s1,s2,s3,s4
ENSG00000000003,100,120,5,8
ENSG00000000419,0,2,300,280
```

```bash
Rscript run_pipeline.R --counts counts.csv --metadata sample_sheet.csv --config config/config.yaml
```

### 4.2 featureCounts (`featurecounts`)

Subread `featureCounts` output: a leading `#` comment line, then
`Geneid Chr Start End Strand Length` and one column per sample (often BAM paths,
which are reduced to sample names). The `Length` column is kept aside for
optional TPM but is not used for DE.

```
# Program:featureCounts v2.0.1; Command:...
Geneid	Chr	Start	End	Strand	Length	s1.bam	s2.bam
ENSG00000000003	chr1	100	150	+	1001	100	120
```

### 4.3 STAR (`star`)

A **directory** of per-sample `*ReadsPerGene.out.tab` files. STAR prepends four
summary rows (`N_unmapped`, `N_multimapping`, `N_noFeature`, `N_ambiguous`),
which are dropped. Columns 2/3/4 are unstranded / forward / reverse counts;
choose with `--star-strand` (default 2 = unstranded).

```bash
Rscript run_pipeline.R --counts star_dir/ --metadata sample_sheet.csv \
  --config config/config.yaml --star-strand 2
```

### 4.4 HTSeq-count (`htseq`)

A **directory** of per-sample 2-column files (`gene<TAB>count`) with trailing
`__no_feature`, `__ambiguous`, … rows (which are dropped). Sample IDs come from
the file names.

### 4.5 Salmon / kallisto (`salmon`)

A **directory** of per-sample subdirectories each containing `quant.sf`
(Salmon) or `abundance.tsv` (kallisto). Requires a transcript→gene map via
`--tx2gene` (a 2-column CSV: transcript_id, gene_id). Transcript-level
quantifications are summarised to gene level with `tximport`
(`countsFromAbundance="no"`) and the fractional estimated counts are rounded to
integers.

```bash
Rscript run_pipeline.R --counts salmon_dir/ --metadata sample_sheet.csv \
  --config config/config.yaml --tx2gene tx2gene.csv
```

### 4.6 (Ranged)SummarizedExperiment (`rangedSE`)

An `.rds` holding a `SummarizedExperiment` / `RangedSummarizedExperiment` (e.g.
`airway`). Counts come from `assay(se)` and — if you omit `--metadata` — the
sample sheet is taken from `colData(se)`.

```bash
Rscript run_pipeline.R --counts experiment.rds --config config/config.yaml
```

---

## 5. The sample sheet

A CSV with one row per sample. Required:

- A **sample-ID** column (default name `sample`; override with `--sample-col`)
  whose values match the count matrix column names **exactly**. A mismatch is a
  hard FAIL that prints the offending IDs in both directions.
- Every variable named in your `design` formula (e.g. `condition`, `batch`).

Example:

```
sample,condition,batch
s1,control,A
s2,control,B
s3,treated,A
s4,treated,B
```

- **Factors & reference level.** Design variables are coerced to factors. The
  contrast's **denominator** is set as the reference level, so
  `log2FoldChange = log2(numerator / denominator)`.
- **Design ↔ config.** `design` and `contrast` live in the config (below). The
  tested variable must be the **last** term of the design.

---

## 6. Configuration reference

All parameters live in `config/config.yaml` (copy from
`config/config.example.yaml`). Every value is echoed to the log at startup, and a
missing required field is a clear startup error. Keep the file ASCII.

| Key | Meaning | Default |
|---|---|---|
| `seed` | Global RNG seed (set before anything stochastic). | `42` |
| `organism` | `human` → org.Hs.eg.db, `mouse` → org.Mm.eg.db. | `human` |
| `input_format` | `auto` or a forced format (§4). | `auto` |
| `design` | R formula; **tested variable last** (e.g. `~ batch + condition`). | — |
| `contrast` | `[factor, numerator, denominator]`; denominator = reference. | — |
| `qc.min_library_size` | WARN if a sample's total counts fall below this. | `1e6` |
| `qc.min_library_size_fail` | FAIL (hard floor) below this. | `1e5` |
| `qc.min_genes_detected` | WARN if a sample expresses fewer genes. | `8000` |
| `qc.min_replicates_per_group` | FAIL if any tested-factor group has fewer. | `2` |
| `qc.min_count` | Low-count filter: keep genes with ≥ this in ≥ min replicates. | `10` |
| `qc.max_pct_genes_filtered` | WARN if the filter removes more than this %. | `90` |
| `qc.cooks_outlier` | Flag Cook's-distance outlier samples (WARN only). | `true` |
| `de.padj_cutoff` | BH-adjusted p significance cut (also DESeq2 `alpha`). | `0.05` |
| `de.lfc_cutoff` | \|log2FC\| threshold for "biologically notable". | `1` |
| `de.shrink` | LFC shrinkage: `apeglm` / `ashr` / `normal` / `none`. | `apeglm` |
| `de.test` | `Wald` (per-coefficient/contrast) or `LRT` (likelihood-ratio, full vs reduced — for multi-level factors / time courses). | `Wald` |
| `de.reduced` | Reduced-model formula string; **required** when `de.test: LRT` (e.g. `"~ batch"`). | `null` |
| `de.contrasts` | Optional list of *additional* results beyond `contrast`: each a `[factor, num, denom]` vector or a `{name, coef}` mapping (a `resultsNames()` coefficient — how you pull an interaction effect). Each gets its own DE table + volcano. | `null` |
| `explore.pca_ntop` | # most-variable genes for PCA. 500 = DESeq2 `plotPCA()` convention (not a hard best practice); 1000–2000 broadens it; a very large value ≈ all genes. | `500` |
| `enrichment.enable` | Run GO ORA + GSEA. | `true` |
| `enrichment.ontology` | GO ontology: `BP` / `MF` / `CC`. | `BP` |
| `enrichment.padj_cutoff` | Significance cut for enriched terms. | `0.05` |
| `outputs.outdir` | Root output directory. | `results/` |

---

## 7. CLI reference

```
Rscript run_pipeline.R [options]
```

| Flag | Description |
|---|---|
| `--counts <file\|dir>` | Path to counts (ignored with `--validate`). |
| `--metadata <csv>` | Sample sheet (optional for `rangedSE`, which carries colData). |
| `--config <yaml>` | Config file (**required**). |
| `--input-format <fmt>` | `auto` (default) or a forced format (§4). |
| `--tx2gene <csv>` | transcript→gene map (**required** for salmon/kallisto). |
| `--outdir <dir>` | Output root (default: `outputs.outdir` from config). |
| `--sample-col <name>` | Sample-ID column in the metadata (default `sample`). |
| `--star-strand <n>` | STAR count column: 2=unstranded, 3=fwd, 4=rev (default 2). |
| `--force` | Downgrade QC FAIL→WARN and continue (logged; stamped on report). |
| `--validate` | Run the airway positive-control self-test (ignores `--counts`). |

**Exit codes:** `0` success; `1` on any QC FAIL (without `--force`) or a
`--validate` assertion failure.

---

## 8. Outputs

Everything is written under `--outdir` (default `results/`):

```
results/
├── report.html                     # self-contained interactive report
├── qc/
│   ├── qc_results.tsv               # every QC check: status, value, message
│   └── variance_vs_metadata.tsv     # R^2 of each PC vs each metadata variable
├── plots/
│   ├── pca.png                      # PCA (static; interactive version in report)
│   ├── pca_batch_corrected.png      # PCA after removeBatchEffect (viz only)
│   ├── variance_vs_metadata.png     # PC-vs-metadata R^2 heatmap
│   ├── volcano.png                  # volcano (static)
│   └── enrichment_dotplot.png       # top enriched GO terms
├── tables/
│   ├── de_results.tsv               # FULL DE table (see columns below)
│   ├── pca_scores.tsv               # PCA coordinates (feeds interactive report)
│   ├── enrichment_ora.tsv           # over-representation analysis
│   └── enrichment_gsea.tsv          # gene-set enrichment analysis
└── logs/
    └── run_<timestamp>.log          # timestamped log + resolved config + sessionInfo()
```

**DE table columns (`de_results.tsv`):**

| Column | Meaning |
|---|---|
| `gene_id` | Input gene identifier (e.g. Ensembl). |
| `symbol` | Gene symbol (mapped via OrgDb; NA if unmapped). |
| `entrez` | Entrez ID (used for enrichment). |
| `baseMean` | Mean of normalized counts across all samples. |
| `log2FoldChange` | Shrunken log2(numerator/denominator). Positive = up in numerator. |
| `lfcSE` | Standard error of the log2 fold change. |
| `stat` | Wald statistic (may be NA after apeglm/ashr shrinkage). |
| `pvalue` | Raw Wald test p-value. |
| `padj` | **Benjamini-Hochberg adjusted p — use this for significance.** |
| `significant` | `TRUE` if `padj < de.padj_cutoff` **and** `\|log2FC\| ≥ de.lfc_cutoff`. |

Genes are ordered by `padj` (most significant first).

---

## 9. Reading the QC report

Open `results/report.html` in any browser (it is fully self-contained — safe to
email or archive). Sections:

- **QC panel** — a red/amber/green table.
  - **Green = PASS**: the check is satisfied.
  - **Amber = WARN**: proceed, but review (e.g. a low-depth sample, a large
    fraction of genes filtered, a candidate outlier). WARNs never stop the run.
  - **Red = FAIL**: a fundamental problem that stops the run unless `--force`
    (misaligned samples, too few replicates, a confounded design, normalized
    input). Hover a check name for what it means.
- **PCA** — samples should separate by the biological condition. Separation
  driven by a technical variable signals a batch effect.
- **Variance explained by metadata** — R² of each PC against each *informative*
  metadata variable (per-sample identifiers are excluded automatically). A
  technical variable dominating PC1 is a batch-effect flag.
- **Differential expression** — summary tiles, an interactive volcano (top genes
  labelled; hover any point; every gene searchable in the table), and the full
  sortable/searchable/exportable DE table.
- **Functional enrichment** — GO ORA dotplot + tables (and GSEA if significant).
- **Reproducibility** — seed, timestamp, organism, batch variables, and where to
  find the log and `renv.lock`.

The console prints the same PASS/WARN/FAIL in color, ending with a boxed summary
banner listing every WARN and FAIL.

---

## 10. Reproducibility model

- **Seed.** `set.seed(config$seed)` runs before any stochastic step.
- **Pinned dependencies.** `renv.lock` records exact versions; `renv::restore()`
  reproduces the environment. `environment.yml` is the conda alternative.
- **Captured environment.** Each run appends `sessionInfo()` to
  `results/logs/run_<timestamp>.log`, and echoes the fully resolved config at the
  top of that log.
- **Config-driven.** No thresholds are hard-coded; a run is fully described by
  its config + lockfile.

**To reproduce a past run:** check out the same commit, `renv::restore()` from the
committed `renv.lock`, and re-run with the same config. The airway `--validate`
pass reproduces identical DEG calls across runs on the same environment.

---

## 10a. Resumability (checkpoint cache)

Every run is **resumable**. The expensive, deterministic stages — the DESeq fit
(size factors + dispersions + GLM), the VST, per-contrast LFC shrinkage +
annotation, and enrichment (ORA/GSEA) — are memoized to `<outdir>/cache/` as
`.rds` files. If a run is interrupted (a crash, an HPC pre-emption, `Ctrl-C`) or
simply re-run, finished stages reload instantly and only what changed is
recomputed.

**How staleness is detected.** Each cache entry's key is a content hash of
everything that determines its result: the count matrix, the sample sheet, the
relevant config values, the upstream stage's key, **and a fingerprint of the
compute code itself** (the deparsed bodies of the compute functions). So:

- Re-running unchanged → every stage reloads; near-instant.
- Change `padj_cutoff` → the DESeq fit and VST are **reused**; only DE +
  enrichment recompute (their keys include the cutoff via `alpha`).
- Edit a compute function (committed *or not*) → its stage and everything
  downstream recompute. Correctness is never traded for reuse.

**Control.**

- On by default. Disable in config with `resume.enable: false`, or per-run with
  the `--no-resume` flag (computes every stage fresh, writes/reads no cache).
- The cache lives under the output directory, so deleting `results/` (or just
  `results/cache/`) clears it. It is git-ignored.

```bash
Rscript run_pipeline.R --counts counts.tsv --metadata samples.csv --config config/config.yaml   # resumable
Rscript run_pipeline.R ... --no-resume                                                            # force fresh
```

On HPC, this means a job killed by a wall-clock limit can simply be resubmitted:
it picks up right after the last completed stage.

---

## 11. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `Sample IDs do not match…` | Count column names ≠ sample-sheet IDs. | Align them exactly; check `--sample-col`. The error lists the offending IDs both ways. |
| `…looks like normalized data (TPM/CPM/FPKM)…` | You supplied normalized values. | Provide raw integer counts. For Salmon/kallisto use `--input-format salmon` + `--tx2gene`. |
| `design … is rank-deficient (confounded covariate)` | A covariate is perfectly confounded with the tested factor (e.g. every treated sample is batch B). | Redesign so batches are crossed with condition, or drop the confounded term from `design`. The printed cross-tab shows the culprit. |
| `group '…' has N replicate(s); need >= 2` | Too few replicates in a group. | Add replicates, or lower `qc.min_replicates_per_group` (not recommended — DE needs within-group replication). |
| `<OrgDb> not installed; skipping…` | Missing annotation DB. | `Rscript -e 'BiocManager::install("org.Hs.eg.db")'` (or `org.Mm.eg.db`). Enrichment/symbols are skipped, not fatal. |
| `Could not auto-detect a supported format` | Ambiguous input. | Pass `--input-format` explicitly (§4). |
| Config parses as empty / all fields "missing" | Non-ASCII characters in `config.yaml` under a C locale. | Keep the config ASCII (the loader reads UTF-8, but ASCII is safest). |
| Report not rendered | `rmarkdown`/`pandoc` missing. | Install pandoc; the analysis outputs (tables/plots) are still written regardless. |
| Re-run reused an old result I didn't expect | The checkpoint cache saw unchanged inputs/config/code. | It only reuses when the content hash matches, but to force a fully fresh run pass `--no-resume` (or delete `<outdir>/cache/`). |

---

## 12. Interpreting results responsibly

- **DEG ≠ causation.** Differential expression is an association under your
  design; it does not establish mechanism or direction of causation.
- **`padj`, not `pvalue`.** Always judge significance on the adjusted p-value.
- **Effect size matters.** A tiny fold change can be "significant" in a
  high-powered design; the `de.lfc_cutoff` guards against over-interpreting it.
- **Enrichment caveats.** ORA depends on your significance cut and the chosen
  gene universe; GSEA depends on the ranking metric. Treat enriched terms as
  hypotheses, not conclusions, and beware annotation bias toward well-studied
  genes.
- **Batch handling.** Batch is corrected **in the DE model** (in the design), not
  by feeding batch-corrected values to DESeq2. The batch-corrected PCA is a
  *visualization* only.
- **"Top" ≠ a specific named gene.** Ranking is by statistical evidence (`padj`);
  a famous marker for your system may be highly significant yet not rank #1.
```
