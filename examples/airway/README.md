# Positive-control output — `airway`

These are the **committed results of the built-in positive-control run** on the
Bioconductor [`airway`](https://bioconductor.org/packages/airway/) dataset
(dexamethasone-treated airway smooth-muscle cells, 8 samples, design
`~ cell + dex`, contrast `dex: trt vs untrt`). They were produced by:

```bash
Rscript run_pipeline.R --validate --config config/config.example.yaml
```

The run recovers the known glucocorticoid response — **CRISPLD2, DUSP1, KLF15 are
significant and up in the treated group** — and the validation exits 0.

| File | What it is |
|---|---|
| [`report.html`](report.html) | The full self-contained interactive report (download & open in a browser). |
| [`qc/qc_results.tsv`](qc/qc_results.tsv) | Every QC check with status/value/message (10 PASS). |
| [`qc/run_metrics.tsv`](qc/run_metrics.tsv) | Run provenance: version, runtime, package versions, counts. |
| [`qc/variance_vs_metadata.tsv`](qc/variance_vs_metadata.tsv) | R² of each PC vs each metadata variable. |
| [`tables/de_results.tsv`](tables/de_results.tsv) | Full DE table (17,199 genes; 831 significant). |
| [`tables/enrichment_ora.tsv`](tables/enrichment_ora.tsv) | Over-represented GO:BP terms (481). |
| [`plots/`](plots/) | PCA, batch-corrected PCA, variance heatmap, volcano, enrichment dotplot. |
| [`logs/run.log`](logs/run.log) | Timestamped run log + resolved config + `sessionInfo()`. |

> `report.html` embeds plotly.js and DataTables, so GitHub shows it as raw HTML
> rather than rendering it inline — download the file (or use a raw-HTML viewer)
> to see the interactive report.
