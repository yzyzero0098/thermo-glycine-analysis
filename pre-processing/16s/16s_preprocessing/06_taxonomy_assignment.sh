#!/bin/bash

echo "Step 06: Taxonomy classification"

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table feature-frequency-filtered-table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.txt \
  --o-visualization taxa-bar-plots.qzv

echo "Taxonomy assignment finished"
