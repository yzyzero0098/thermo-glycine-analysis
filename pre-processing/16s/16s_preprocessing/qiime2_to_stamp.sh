#!/bin/bash

# 1. 작업 폴더 생성
mkdir -p function

# 2. 필요한 QIIME2 artifact 복사
cp rooted-tree.qza function/
cp table.qza function/
cp taxonomy.qza function/
cp *.txt function/

# 3. 작업 폴더 이동
cd function

# 4. QIIME2 artifact export
qiime tools export --input-path table.qza --output-path .
qiime tools export --input-path taxonomy.qza --output-path .

# 5. STAMP format에 맞게 taxonomy header 수정
sed -i '1s/Feature ID/#OTUID/;1s/Taxon/taxonomy/' taxonomy.tsv

# 6. taxonomy metadata를 feature table에 추가
biom add-metadata \
-i feature-table.biom \
-o feature-table-with-taxonomy.biom \
--observation-metadata-fp taxonomy.tsv \
--sc-separated taxonomy

# 7. BIOM → TSV 변환 (STAMP용)
biom convert \
-i feature-table-with-taxonomy.biom \
-o feature-table.tsv \
--to-tsv \
--header-key taxonomy

echo "STAMP input table 생성 완료"
