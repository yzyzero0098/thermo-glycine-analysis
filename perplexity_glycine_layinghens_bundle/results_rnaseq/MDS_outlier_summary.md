# RNA-seq MDS outlier check

- Method: edgeR `plotMDS` on TMM-normalized counts after `filterByExpr` filtering
- Outlier candidate rule: within-group distance-to-centroid robust z-score > 2.5 on MDS1/MDS2
- Note: this is a screening rule, not an automatic exclusion criterion

## Candidate samples
- P1-2-JM-NC-R (NC, distance=0.986, robust_z=3.60)
- P1-2-JM-PC-R (PC, distance=0.659, robust_z=10.79)
- P1-5-JM-PC-R (PC, distance=0.657, robust_z=10.73)
