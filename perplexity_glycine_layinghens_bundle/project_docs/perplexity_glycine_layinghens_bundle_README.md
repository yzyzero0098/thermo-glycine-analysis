# Perplexity Bundle for `glycine_layinghens`

This bundle is a curated project package for literature-guided interpretation, hypothesis refinement, and additional analysis planning in the `glycine_layinghens` project.

It excludes raw FASTQ/BAM files to keep the package practical for reading tools, but includes the main inputs and result summaries needed for:

- project background review
- 16S microbiome interpretation
- PICRUSt2 functional interpretation
- RNA-seq DEG interpretation
- phenotype caveat assessment
- additional analysis planning

## Directory guide

- `project_docs`
  - proposal/manuscript context and presentation materials
- `inputs_16s`
  - key QIIME2/PICRUSt2 input and output tables
- `inputs_rnaseq`
  - RNA-seq preprocessing summary, count matrix, sample metadata, and reference annotation
- `results_16s`
  - phyloseq, DA, and functional-analysis summaries
- `results_rnaseq`
  - RNA-seq DEG summaries, robustness comparisons, and MDS outlier check

## Notes for interpretation

- `PICRUSt2` outputs are predicted functions, not directly measured metagenomics.
- phenotype integration should be treated cautiously because some phenotype tables may contain repeated group means rather than true sample-level raw measurements.
- RNA-seq DEG counts differ by workflow; both `DESeq2` and `edgeR` outputs are included, along with a legacy `edgeR` reanalysis for comparison.

## Suggested priority files

If a reading tool cannot process everything, start with:

1. `project_docs/glycine_길교수님방_publish.pdf`
2. `results_16s/glycine_16S_interpretation_report.md`
3. `results_16s/analysis_report.md`
4. `results_16s/DA_cross_method_summary.md`
5. `results_16s/functional_analysis_verified_summary.md`
6. `results_rnaseq/RNAseq_functional_analysis_report.md`
7. `results_rnaseq/DEG_method_comparison.csv`
8. `results_rnaseq/MDS_outlier_summary.md`
9. `inputs_16s/path_abun_unstrat.tsv.gz`
10. `inputs_rnaseq/merged_count_matrix_with_annotation.csv`
