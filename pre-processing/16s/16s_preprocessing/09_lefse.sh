#!/bin/bash

echo "Running LEfSe"

mkdir -p lefse
cd lefse

format_input.py feature-table.tsv lefse.in \
 -c 1 \
 -u 2 \
 -o 1000000


run_lefse.py lefse.in lefse.res


plot_res.py lefse.res lefse_LDA.png \
 --format png


plot_cladogram.py lefse.res cladogram.png

echo "LEfSe done"
