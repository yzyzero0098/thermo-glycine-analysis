# Reproducibility Notes

This repository is a curated, lightweight reproduction package for the thermo-glycine analysis. It contains analysis scripts, interpretation notes, and selected compact result tables. It does not contain raw sequencing reads or large intermediate files.

## Intended Reproduction Level

The repository supports three levels of reuse:

1. Review the analysis logic from scripts and notes.
2. Reproduce downstream figures/tables when matching derived input tables are available.
3. Re-run the full workflow if raw data are obtained separately and placed into the expected local paths.

## Data Not Included

The following files are intentionally excluded:

- Raw FASTQ files
- BAM/BAI/SAM alignment outputs
- QIIME2 artifacts and visualizations (`.qza`, `.qzv`)
- Large compressed archives
- PowerPoint, Word, PDF, and manuscript working files
- Private or unpublished raw data

If this repository is used for a paper, raw data should be deposited separately in an appropriate repository such as NCBI SRA, ENA, DDBJ, PRIDE, MetaboLights, Figshare, Zenodo, or an institutional repository depending on data type and publication policy.

## Analysis Components

### 16S / microbiome

Main files:

- `pre-processing/16s/16s_preprocessing/*.sh`
- `pre-processing/16s/16s_preprocessing/rarefaction_depth.py`
- `pre-processing/16s/16s_preprocessing/nsti_histogram.R`
- `Functional analysis/Glycine_HS_16S_Pipeline_Complete.R`
- `Functional analysis/analysis_part1_diversity.py`
- `Functional analysis/analysis_part2_taxa_DA.py`

Typical workflow:

1. Import reads and metadata into QIIME2.
2. Denoise with DADA2.
3. Filter feature table and build phylogeny.
4. Run alpha/beta diversity analyses.
5. Assign taxonomy and export tables.
6. Run PICRUSt2 and LEfSe-style downstream summaries where applicable.
7. Recreate diversity, composition, differential abundance, and dose-response summaries.

### RNA-seq / transcriptome

Main files:

- `pre-processing/rna_seq/*.R`
- `pre-processing/rna_seq/1.fastqc/xx-00-qcext.sh`
- `pre-processing/rna_seq/3.trimQC/xx-01-qcext.sh`
- `glycine/tables/rnaseq/*.csv`

Typical workflow:

1. Run raw read QC.
2. Trim reads and inspect post-trim QC.
3. Align or pseudoalign reads against the chosen reference genome and annotation.
4. Generate count matrices.
5. Run differential expression and sensitivity/outlier scans.
6. Summarize MDS/PCA, DEG counts, and pathway/GSEA outputs.

### Phenotype and integration

Main files:

- `phenotype/analysis_part3_phenotype.py`
- `glycine/tables/integration/*.csv`
- `glycine/docs/glycine-procrustes-integration.md`
- `glycine/docs/glycine-omics-reanalysis.md`

Typical workflow:

1. Standardize treatment and sample identifiers.
2. Match phenotype, microbiome, and host transcriptome samples.
3. Compute group summaries and distance-based integration summaries.
4. Interpret treatment relationships with attention to matched sample limitations.

## Environment Notes

Exact environment files were not recovered in this curated copy. Recommended modern environment setup:

- QIIME2 for 16S preprocessing
- R 4.x with tidyverse, vegan, phyloseq, DESeq2, edgeR, limma, fgsea as needed
- Python 3.10+ with pandas, numpy, scipy, scikit-bio or scikit-learn, matplotlib/seaborn as needed
- Nextflow/nf-core may be used for a cleaner future RNA-seq or metagenomics rerun

## Suggested Citation / Provenance Text

This repository is a curated analysis-code and derived-summary package prepared from a local thermo-glycine research workspace. Raw data and large intermediate artifacts were excluded from GitHub and should be archived separately for full computational reproducibility.

