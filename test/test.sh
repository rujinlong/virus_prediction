#!/bin/bash

DATA_RAW=$1
DATA_CLEAN=d02_data_clean
ASB=p01_assembly

rm -rf $DATA_CLEAN $ASB
mkdir $DATA_CLEAN $ASB

cd $DATA_RAW
ls *.fq.gz | cut -d'_' -f1 | sort -u > ../ids
cd ../$DATA_CLEAN
cat ../ids | parallel ../preprocess.sh -c qc -1 ../${DATA_RAW}/{}_R1.fq.gz -2 ../${DATA_RAW}/{}_R2.fq.gz -f {} -t 4

cd ../${ASB}
cat ../ids | parallel ../preprocess.sh -c assemble_megahit -1 ../${DATA_CLEAN}/{}_deduped.fq.gz -s ../${DATA_CLEAN}/{}_dedupedsingle.fq.gz -f {} -t 10 -m 0.4

cd ..
cat ${ASB}/*.fna > contigs.fna
nextflow main.nf --fid test --contigs contigs.fna -resume

if [ -e result/test/p05_taxonomy/taxonomy.tsv ];then
    echo "Sucessful"
else
    echo "Failed"
fi
