#!/bin/bash

mkdir -p function

cp table.qza function/
cp taxonomy.qza function/
cp *.txt function/

cd function

echo "Exporting QIIME artifacts"

qiime tools export \
 --input-path table.qza \
 --output-path .


qiime tools export \
 --input-path taxonomy.qza \
 --output-path .


sed -i '1s/Feature ID/#OTUID/;1s/Taxon/taxonomy/' taxonomy.tsv


biom add-metadata \
 -i feature-table.biom \
 -o feature-table-with-taxonomy.biom \
 --observation-metadata-fp taxonomy.tsv \
 --sc-separated taxonomy


biom convert \
 -i feature-table-with-taxonomy.biom \
 -o feature-table.tsv \
 --to-tsv \
 --header-key taxonomy

cd ..

mkdir picrust2

qiime tools export \
 --input-path rep-seqs.qza \
 --output-path picrust2

qiime tools export \
 --input-path table.qza \
 --output-path picrust2

cd picrust2
mv dna-sequences.fasta rep-seqs.fasta


echo "STAMP table ready"
