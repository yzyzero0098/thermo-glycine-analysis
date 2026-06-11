#!/bin/bash

THREADS=64

echo "Step 2: DADA2"

qiime dada2 denoise-paired \
 --i-demultiplexed-seqs paired-end-demux.qza \
 --p-trim-left-f 5 \
 --p-trim-left-r 5 \
 --p-trunc-len-f 260 \
 --p-trunc-len-r 220 \
 --o-table table.qza \
 --o-representative-sequences rep-seqs.qza \
 --o-denoising-stats denoising-stats.qza \
 --p-n-threads $THREADS


qiime feature-table summarize \
 --i-table table.qza \
 --o-visualization table.qzv \
 --m-sample-metadata-file sample-metadata.txt


qiime feature-table tabulate-seqs \
 --i-data rep-seqs.qza \
 --o-visualization rep-seqs.qzv


qiime metadata tabulate \
 --m-input-file denoising-stats.qza \
 --o-visualization denoising_stats.qzv

echo "DADA2 complete"
