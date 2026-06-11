# Thermo Glycine Analysis

Curated analysis workspace for a heat-stress glycine supplementation project in laying hens.

This repository is organized for sharing analysis code, reusable summaries, and selected result tables. Large raw sequencing files and binary research artifacts are intentionally excluded from version control.

## Who Can Reuse This Repository?

This repository is intended for:

- researchers studying nutritional mitigation of heat stress in poultry
- microbiome and transcriptome analysts building reproducible multi-omics workflows
- students learning how to organize 16S, RNA-seq, phenotype, and integration scripts
- maintainers who need a compact public example of omics analysis curation
- reviewers who want to inspect the analysis logic without downloading raw data

The repository is not a raw-data archive. It is a reproducibility layer: code,
lightweight derived tables, workflow notes, and interpretation memos that can be
adapted to related animal, nutrition, microbiome, or host-response studies.

## Open Research Value

This project is shared as an open research-automation and reproducibility
workspace. The main goal is to make the analysis structure reusable beyond the
original study by showing:

- how 16S microbiome, RNA-seq, and phenotype outputs can be organized together
- how treatment contrasts and dose-response logic can be documented
- how lightweight derived outputs can be shared while keeping raw data separate
- how manuscript-facing interpretation notes can remain linked to executable code

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

## Installation And Reproducibility

Clone the repository:

```powershell
git clone https://github.com/yzyzero0098/thermo-glycine-analysis.git
cd thermo-glycine-analysis
```

Recommended local tools:

- R 4.3 or newer for ecological statistics and visualization scripts
- Python 3.10 or newer for lightweight table utilities
- QIIME2 for 16S preprocessing commands
- PICRUSt2 and LEfSe when reproducing functional prediction steps
- standard RNA-seq command-line tools for QC and count preparation

Suggested reproduction order:

1. Read `glycine/README.md` for the curated project map.
2. Review `REPRODUCIBILITY.md` for excluded data and expected inputs.
3. Run preprocessing scripts only after placing raw inputs in a private data directory.
4. Run downstream R scripts against lightweight derived tables where available.
5. Compare generated outputs with the selected summary tables in `glycine/tables/`.

Paths in older scripts may reflect the original local workstation. Treat those
as examples and update them to your own private input/output folders.

## Citation

If this workspace helps your analysis, cite the repository using
`CITATION.cff` or the GitHub citation button.

## Contributing

Issues and pull requests are welcome for reproducibility fixes, documentation
improvements, and portable workflow examples. See `CONTRIBUTING.md`.

## Notes

This is a curated copy from the local research workspace, prepared for GitHub publication. It does not contain the full raw analysis directory.
