#!/bin/bash

INPUT=/PATH/TO/samples_list.csv
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read -r SRA SAMPLE ACCESSION;
do
bash ./backmapper.sh $SAMPLE $ACCESSION $SRA </dev/null
done < $INPUT
IFS=$OLDIFS
