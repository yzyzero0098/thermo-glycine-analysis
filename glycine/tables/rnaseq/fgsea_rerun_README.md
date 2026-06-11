# Re-run fgsea outputs and plots

This set was re-generated in the lab-meeting asset folder using R/fgsea.

## Ranking input
- `NC_vs_G25`, `NC_vs_G50`, `NC_vs_PC`: archived DESeq2 full result tables
- `G50_vs_PC`: saved PyDESeq2 direct-contrast full table
- Ranking statistic: Wald `stat`

## Gene sets
- `msigdbr` Gallus gallus C5 GO:BP
- gene-set size filter: 10 to 500 genes

## Generated files
- full fgsea result tables
- top30 significant tables
- per-contrast top-pathway barplots
- per-contrast top positive/negative enrichment plots
- cross-contrast summary dotplot
- interpretation summary table

## Script
- `experiments/glycine-gsea-rerun/run_fgsea_with_plots.R`
