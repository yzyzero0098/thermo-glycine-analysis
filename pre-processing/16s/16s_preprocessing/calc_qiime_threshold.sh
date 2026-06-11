#!/bin/bash

QZV=$1

if [ -z "$QZV" ]; then
  echo "Usage: ./calc_qiime_threshold.sh table.qzv"
  exit 1
fi

echo "===================================="
echo "Scanning QZV structure"
echo "===================================="

# 내부 파일 자동 탐색
SAMPLE_FILE=$(unzip -l "$QZV" | awk '{print $4}' | grep -i sample | grep -i freq | head -n1)
FEATURE_FILE=$(unzip -l "$QZV" | awk '{print $4}' | grep -i feature | grep -i freq | head -n1)

if [ -z "$SAMPLE_FILE" ]; then
  echo "ERROR: sample frequency file not found in qzv"
  exit 1
fi

if [ -z "$FEATURE_FILE" ]; then
  echo "ERROR: feature frequency file not found in qzv"
  exit 1
fi

echo "Sample file detected : $SAMPLE_FILE"
echo "Feature file detected: $FEATURE_FILE"

echo ""
echo "===================================="
echo "Extracting data"
echo "===================================="

unzip -p "$QZV" "$SAMPLE_FILE" > sample_depth.csv
unzip -p "$QZV" "$FEATURE_FILE" > feature_freq.csv

echo ""
echo "===================================="
echo "Sample depth statistics"
echo "===================================="

awk -F',' '
NR>1{
sum+=$2
if(min=="" || $2<min) min=$2
if($2>max) max=$2
count++
}
END{
print "Samples:",count
print "Min depth:",min
print "Max depth:",max
if(count>0){print "Mean depth:",sum/count}
}' sample_depth.csv

echo ""
echo "Sorting depth values..."

cut -d',' -f2 sample_depth.csv | tail -n +2 | sort -n > depth_sorted.txt

N=$(wc -l depth_sorted.txt | awk '{print $1}')

if [ "$N" -eq 0 ]; then
  echo "ERROR: No sample depth values found"
  exit 1
fi

IDX=$(awk -v n=$N 'BEGIN{printf("%d", n*0.1)}')
if [ "$IDX" -lt 1 ]; then IDX=1; fi

RARE=$(sed -n "${IDX}p" depth_sorted.txt)

echo ""
echo "===================================="
echo "Recommended rarefaction depth"
echo "===================================="
echo "10 percentile depth: $RARE"

echo ""
echo "===================================="
echo "Feature frequency statistics"
echo "===================================="

cut -d',' -f2 feature_freq.csv | tail -n +2 | sort -n > feature_sorted.txt

N=$(wc -l feature_sorted.txt | awk '{print $1}')

IDX=$(awk -v n=$N 'BEGIN{printf("%d", n*0.01)}')
if [ "$IDX" -lt 1 ]; then IDX=1; fi

FREQ=$(sed -n "${IDX}p" feature_sorted.txt)

echo ""
echo "Recommended feature filtering"
echo "min-frequency: $FREQ"
echo "min-samples: 3"

echo ""
echo "===================================="
echo "Recommended QIIME2 pipeline"
echo "===================================="

echo ""
echo "qiime feature-table filter-features \\"
echo "  --i-table table.qza \\"
echo "  --p-min-frequency $FREQ \\"
echo "  --p-min-samples 3 \\"
echo "  --o-filtered-table filtered-table.qza"

echo ""
echo "qiime diversity core-metrics \\"
echo "  --i-table filtered-table.qza \\"
echo "  --p-sampling-depth $RARE \\"
echo "  --m-metadata-file sample-metadata.txt \\"
echo "  --output-dir core-metrics"