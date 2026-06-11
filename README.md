# Thermo Glycine Analysis

Curated analysis workspace for a heat-stress glycine supplementation project in laying hens.

This repository is organized for sharing analysis code, reusable summaries, and selected result tables. Large raw sequencing files and binary research artifacts are intentionally excluded from version control.

## Project Scope

The project combines:

- 16S rRNA microbiome preprocessing and downstream ecological analysis
- RNA-seq quality control, differential expression, and pathway analysis
- Phenotype and host-microbiome integration summaries
- Glycine dose-response and PC-reference comparison notes

Treatment labels used across the analysis:

- `NC`: heat-stress baseline condition
- `PC`: positive-control / energy-substitution reference condition
- `G25`: low-dose glycine supplementation
- `G50`: high-dose glycine supplementation

## Repository Layout

```text
Functional analysis/
  Python and R scripts for 16S diversity, taxa, and dose-response analyses.

phenotype/
  Phenotype integration script.

pre-processing/
  QIIME2, PICRUSt2, LEfSe, RNA-seq QC, and RNA-seq exploratory scripts.

glycine/docs/
  Interpretation notes, manuscript framing, GSEA notes, and integration memos.

glycine/tables/
  Selected compact result tables for 16S, RNA-seq, and integration analyses.

glycine/additional_run/
  Follow-up PICRUSt and targeted integration notes and summary tables.
```

## Data Availability

Large and raw files are not included here. Excluded materials include:

- FASTQ, BAM, BAI, QZA, QZV, BIOM, and compressed raw outputs
- MZmine installers and other software binaries
- large zip archives
- PowerPoint, Word, PDF, and other presentation/manuscript working files
- intermediate HTML, log, and generated visualization folders

The repository is intended to document analysis logic and share lightweight derived outputs. Raw data should be stored separately according to lab, institution, or journal data-sharing requirements.

## Suggested Starting Points

- `glycine/README.md`
- `glycine/docs/glycine-omics-reanalysis.md`
- `glycine/docs/glycine-procrustes-integration.md`
- `Functional analysis/Glycine_HS_16S_Pipeline_Complete.R`
- `pre-processing/16s/16s_preprocessing/01_qiime_import.sh`
- `pre-processing/rna_seq/rna_seq_visualization_glycine.R`

## Notes

This is a curated copy from the local research workspace, prepared for GitHub publication. It does not contain the full raw analysis directory.
