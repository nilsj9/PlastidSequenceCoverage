#!/bin/bash

NCBIENT=/PATH/TO/NCBI-ENTREZ-DIRECT
PACVR=/PATH/TO/PACVRRSCRIPT
SAMPLELOC=/PATH/TO/SAMPLES

# DEFINITIONS
SAMPLE=$1
ACCESSION=$2
SRA=$3

# conduct PACVr Rscript
Rscript $PACVR -k $SAMPLELOC/$SAMPLE/${ACCESSION}.gb -b $SAMPLELOC/$SAMPLE/${ACCESSION}_mapping_OneMoreLocations.sorted.bam -t 0.5 -r TRUE -v TRUE -o $SAMPLELOC/$SAMPLE/${ACCESSION}_CoverageDepth.pdf

# merge metadata to one csv
${NCBIENT}esearch -db sra -query $SRA | efetch -format runinfo | cut -f20,7 -s -d, > $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/${ACCESSION}_SRA_meta.csv 
${NCBIENT}esearch -db nuccore -query $ACCESSION | efetch -format gb | sed -n '/Assembly Method/,/Sequencing/{/Sequencing/!p;}'  | awk '{$1=$1;print}' | sed 's/:: /\n/g' | paste -s -d ' ' | sed 's/  /\n/g' | awk /./ > $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/${ACCESSION}_assembly_tech.csv
find $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/ -type f -size 0 -exec rm {} \;
paste -d, $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/*.csv > $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/${ACCESSION}_metadata.csv