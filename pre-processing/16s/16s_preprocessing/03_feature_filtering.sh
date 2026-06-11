#!/bin/bash

echo "Step 03: Feature filtering"

qiime feature-table filter-samples \
  --i-table table.qza \
  --p-min-frequency 1500 \
  --o-filtered-table sample-frequency-filtered-table.qza

qiime feature-table filter-features \
  --i-table sample-frequency-filtered-table.qza \
  --p-min-frequency 31 \
  --p-min-samples 3 \
  --o-filtered-table feature-frequency-filtered-table.qza

qiime feature-table summarize \
  --i-table feature-frequency-filtered-table.qza \
  --o-visualization feature-frequency-filtered-table.qzv \
  --m-sample-metadata-file sample-metadata.txt

echo "Feature filtering done"
