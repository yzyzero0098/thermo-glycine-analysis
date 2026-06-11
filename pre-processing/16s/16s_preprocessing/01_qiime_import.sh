#!/bin/bash

echo "Step 1: Import FASTQ"

qiime tools import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path importing_file.csv \
 --output-path paired-end-demux.qza \
 --input-format PairedEndFastqManifestPhred33


qiime demux summarize \
 --i-data paired-end-demux.qza \
 --o-visualization demux.qzv

echo "Import complete"
