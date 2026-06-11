#!/bin/bash

SAMPLING_DEPTH=4030

echo "Step 05: Diversity analysis"

qiime diversity alpha-rarefaction \
  --i-table feature-frequency-filtered-table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth $SAMPLING_DEPTH \
  --m-metadata-file sample-metadata.txt \
  --o-visualization alpha-rarefaction_${SAMPLING_DEPTH}.qzv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table feature-frequency-filtered-table.qza \
  --p-sampling-depth $SAMPLING_DEPTH \
  --m-metadata-file sample-metadata.txt \
  --output-dir core-metrics-results

echo "Diversity analysis finished"
