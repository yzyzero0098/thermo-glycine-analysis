#!/bin/bash

BasePath=$(pwd)
#need to check PATH!!!!!!!!!!!!!!!!!!
qcPath=${BasePath}/

cd ${qcPath}

#unzip
for i in $(ls *.zip)
do
  unzip ${i}
  echo ${i}
done

#extraction of raw sequence & gc percentage
[[ -f ${qcPath}/1.SAMPLENAME ]] && rm -f ${qcPath}/1.SAMPLENAME
[[ ! -f ${qcPath}/1.SAMPLENAME ]] && touch ${qcPath}/1.SAMPLENAME

[[ -f ${qcPath}/2.SEQUENCE ]] && rm -f ${qcPath}/2.SEQUENCE
[[ ! -f ${qcPath}/2.SEQUENCE ]] && touch ${qcPath}/2.SEQUENCE

[[ -f ${qcPath}/3.QCPERCENT ]] && rm -f ${qcPath}/3.QCPERCENT
[[ ! -f ${qcPath}/3.QCPERCENT ]] && touch ${qcPath}/3.QCPERCENT

for f in $(ls -l | grep '^d' | awk '{print $NF}')
do
  cd ${f}
  echo ${f}
  awk 'NR==4 {print $2}' fastqc_data.txt >> ${qcPath}/1.SAMPLENAME
  awk 'NR==7 {print $3}' fastqc_data.txt >> ${qcPath}/2.SEQUENCE
  awk 'NR==10 {print $2}' fastqc_data.txt >> ${qcPath}/3.QCPERCENT
  cd ..
done

paste 1.SAMPLENAME 2.SEQUENCE 3.QCPERCENT > 4.FINALEXTRACTION
