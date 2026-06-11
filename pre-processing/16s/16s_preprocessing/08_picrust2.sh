#!/bin/bash

echo "=============================="
echo "Running PICRUSt2"
echo "=============================="

cd picrust2

# run PICRUSt2
picrust2_pipeline.py \
 -s rep-seqs.fasta \
 -i feature-table.biom \
 -o picrust2_out \
 -p 32

echo ""
echo "PICRUSt2 finished"
echo "=============================="

############################################
# QC summary
############################################

echo ""
echo "PICRUSt2 QC summary"
echo "------------------------------"

# total ASV
TOTAL_ASV=$(grep -c ">" rep-seqs.fasta)
echo "Total ASVs input: $TOTAL_ASV"

# predicted KO table 존재 확인
if [ -f picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz ]; then
    echo "KO prediction table generated"
else
    echo "WARNING: KO table missing"
fi

# pathway table 존재 확인
if [ -f picrust2_out/pathways_out/path_abun_unstrat.tsv.gz ]; then
    echo "Pathway prediction table generated"
else
    echo "WARNING: pathway table missing"
fi

# NSTI statistics
if [ -f picrust2_out/marker_predicted_and_nsti.tsv.gz ]; then
    echo ""
    echo "NSTI statistics:"
    zcat picrust2_out/marker_predicted_and_nsti.tsv.gz | \
    awk 'NR>1 {sum+=$2; n++} END {print "Mean NSTI:", sum/n}'
fi

# KO table의 총 줄 수 (functional feature의 개수)
echo ""
echo "Functional feature summary"

KO_COUNT=$(zcat picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz | wc -l)

echo "Total KO features: $KO_COUNT"

# Pathway의 개수
PATHWAY_COUNT=$(zcat picrust2_out/pathways_out/path_abun_unstrat.tsv.gz | tail -n +2 | wc -l)
echo "Total pathways predicted: $PATHWAY_COUNT"

echo ""
echo "Generating NSTI histogram"

Rscript ../nsti_histogram.R

echo ""
echo "PICRUSt2 pipeline completed successfully"
echo "=============================="

