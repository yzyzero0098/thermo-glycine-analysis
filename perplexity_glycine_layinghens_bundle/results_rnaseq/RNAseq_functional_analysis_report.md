# RNA-seq Functional Analysis

- Input counts: featureCounts outputs from Gallus 115 remapping
- Samples: 25
- Genes after low-count filtering: 17421

## DEG summary
- NC_vs_PC: total 86, up 59, down 27
- NC_vs_G25: total 4, up 1, down 3
- NC_vs_G50: total 884, up 818, down 66

## DEG method comparison
- NC_vs_PC / DESeq2: tested 17396, total 86, up 59, down 27
- NC_vs_PC / edgeR: tested 17421, total 44, up 26, down 18
- NC_vs_G25 / DESeq2: tested 17396, total 4, up 1, down 3
- NC_vs_G25 / edgeR: tested 17421, total 0, up 0, down 0
- NC_vs_G50 / DESeq2: tested 17396, total 884, up 818, down 66
- NC_vs_G50 / edgeR: tested 17421, total 944, up 889, down 55

## Treatment colors
- NC: #1976D2
- PC: #388E3C
- G25: #F57C00
- G50: #D32F2F

## Notes
- DEG significance threshold: FDR < 0.05 and |log2FC| >= 1.
- DESeq2 results were exported with independentFiltering = FALSE for fair cross-comparison of DEG counts.
- edgeR analysis used TMM normalization and quasi-likelihood GLM testing.
- GO enrichment used Biological Process terms.
- KEGG enrichment was attempted through clusterProfiler/KEGGREST and may depend on remote KEGG access.
