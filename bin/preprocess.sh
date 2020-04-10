#!/bin/bash

DIR_VIROMEQC=${HOME}/usr/opt/viromeqc

while getopts "1:2:c:f:m:s:t:w:" arg
do
    case $arg in
        c)
          cmd=$OPTARG
          ;;
        1)
          R1=$OPTARG
          Rint=$OPTARG
          ;;
        2)
          R2=$OPTARG
          ;;
        f)
          fid=$OPTARG
          ;;
        m)
          mem=$OPTARG
          ;;
        s)
          Rs=$OPTARG
          ;;
        t)
          t=$OPTARG
          ;;
        w)
          enrichment=$OPTARG
          ;;
        ?)
          echo "unknow argument"
          exit 1
          ;;
    esac
done


if [ $cmd = 'qc' ];then
    # -1=$R1, -2=$R2, -f=$fid, -t=$t
    fastp -i $R1 -I $R2 -o ${fid}_clean_R1.fq.gz -O ${fid}_clean_R2.fq.gz --unpaired1 ${fid}_singletons.fq.gz --unpaired2 ${fid}_singletons.fq.gz --failed_out ${fid}_fail.fq.gz -p -w $t -n 1 -l 30 -5 -W 4 -M 20 -r -c -g -j ${fid}.json -h ${fid}.html
    dedupe.sh in1=${fid}_clean_R1.fq.gz in2=${fid}_clean_R2.fq.gz out=${fid}_deduped.fq.gz outd=${fid}_dup.fq.gz ac=f threads=$t
    dedupe.sh in=${fid}_singletons.fq.gz out=${fid}_dedupedsingle.fq.gz outd=${fid}_dupsingle.fq.gz ac=t threads=$t
elif [ $cmd = 'viromeqc' ];then
    # -1=$Rint, -s=$Rs, -f=$fid, -w=$enrichment {human,environmental}
    ${DIR_VIROMEQC}/viromeQC.py -i $Rint $Rs -o vqc_report_${fid}.tsv -w $enrichment
elif [ $cmd = 'assemble_spades' ];then
    # -1=$Rint, -s=$Rs, -f=$fid, -t=$t, -m=$mem
    metaspades.py -o asb_${fid} --12 $Rint -s $Rs -t $t -m $mem
    sed "s/^>/>${fid}_/" asb_${fid}/contigs.fasta | sed 's/ /_/g' > ${fid}.fna
elif [ $cmd = 'assemble_megahit' ];then
    # -1=$Rint, -s=$Rs, -f=$fid, -t=$t, -m=$mem
    megahit -o asb_${fid} --12 $Rint -r $Rs -t $t -m $mem
    sed "s/^>/>${fid}_/" asb_${fid}/final.contigs.fa | sed 's/ /_/g' > ${fid}.fna
fi
