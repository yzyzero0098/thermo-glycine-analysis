#!/bin/bash

echo "Network analysis"

mkdir -p network
cd network

# SparCC correlation

sparcc.py feature-table.tsv \
 --cor_file=sparcc_cor.tsv \
 --cov_file=sparcc_cov.tsv


echo "Network file ready for Cytoscape"
