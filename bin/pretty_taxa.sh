#!/bin/bash

PREDICTED_TAXA=$1
NAME2TAXAID=nf_tmp_name2taxid.tsv
TAXAID2LINEAGE=nf_tmp_taxaid2lineage.tsv

rev $PREDICTED_TAXA | cut -f1 | rev | sed '1d' | sort -u | taxonkit name2taxid > $NAME2TAXAID
sleep 2
cut -f2 $NAME2TAXAID | sort -u | sed '/^$/d' | taxonkit lineage | awk '$2!=""' | taxonkit reformat | cut -f 1,3 | sed -r 's/;/\t/g' | csvtk add-header -t -n taxid,Superkingdom,Phylum,Class,Order,Family,Genus,Species > $TAXAID2LINEAGE
sleep 2
 
if [ -e "$TAXAID2LINEAGE" ];then
    echo "Success"
else
    echo "Something wrong!"
fi

# rm $NAME2TAXAID $TAXAID2LINEAGE
