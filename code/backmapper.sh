#!/bin/bash

# REQUIREMENTS
#ncbi-entrez-direct 
#sra-toolkit >2.10.8 
#samtools >1.10
#BOWTIE2 >2.3.4.1
#TRIMMOMATIC
#RCircos
#PACVr
#genbankr

# SET PATHS
SAM=/PATH/TO/SAMTOOLS
SRAT=/PATH/TO/SRA-TOOLKIT
NCBIENT=/PATH/TO/NCBI-ENTREZ-DIRECT
BOWTIE=/PATH/TO/BOWTIE2
TRIMMOMATIC=/PATH/TO/TRIMMOMATIC
PACVR=/PATH/TO/PACVRRSCRIPT
SAMPLELOC=/PATH/TO/SAMPLES

# DEFINITIONS
SAMPLE=$1
ACCESSION=$2
SRA=$3

# FOLDERS
mkdir -p $SAMPLELOC/$SAMPLE/Backmapping

# download GenBank
${NCBIENT}esearch -db nucleotide -query $ACCESSION | ${NCBIENT}efetch -format gb > $SAMPLELOC/$SAMPLE/${ACCESSION}.gb

# download reference fasta
${NCBIENT}esearch -db nucleotide -query $ACCESSION | ${NCBIENT}efetch -format fasta > $SAMPLELOC/$SAMPLE/${ACCESSION}_ref.fasta

# download reads
cd $SAMPLELOC/$SAMPLE
${SRAT}prefetch --max-size 50000000 $SRA
cd /scratch/nilsj24/
${SRAT}fasterq-dump.2.10.8 --split-3 --skip-technical $SAMPLELOC/$SAMPLE/$SRA/$SRA.sra -O $SAMPLELOC/$SAMPLE/Backmapping

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE $SAMPLELOC/$SAMPLE/Backmapping/${SRA}.sra_1.fastq $SAMPLELOC/$SAMPLE/Backmapping/${SRA}.sra_2.fastq  $SAMPLELOC/$SAMPLE/Backmapping/${SRA}_1_PE.fastq $SAMPLELOC/$SAMPLE/Backmapping/${SRA}_1_SE.fastq $SAMPLELOC/$SAMPLE/Backmapping/${SRA}_2_PE.fastq $SAMPLELOC/$SAMPLE/Backmapping/${SRA}_2_SE.fastq ILLUMINACLIP:/$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

# conduct backmapping
REF=$SAMPLELOC/$SAMPLE/${ACCESSION}_ref.fasta
READSPER1=$SAMPLELOC/$SAMPLE/Backmapping/${SRA}_1_PE.fastq
READSPER2=$SAMPLELOC/$SAMPLE/Backmapping/${SRA}_2_PE.fastq

mkdir -p $SAMPLELOC/$SAMPLE/Backmapping/db
${BOWTIE}bowtie2-build $REF $SAMPLELOC/$SAMPLE/Backmapping/db/$ACCESSION
${BOWTIE}bowtie2 -x $SAMPLELOC/$SAMPLE/Backmapping/db/$ACCESSION -1  $READSPER1 -2 $READSPER2 -S $SAMPLELOC/$SAMPLE/Backmapping/${ACCESSION}_mapping.sam
${SAM}samtools view -Sb -F 0x04 $SAMPLELOC/$SAMPLE/Backmapping/${ACCESSION}_mapping.sam > $SAMPLELOC/$SAMPLE/Backmapping/${ACCESSION}_mapping_OneMoreLocations.bam
${SAM}samtools sort $SAMPLELOC/$SAMPLE/Backmapping/${ACCESSION}_mapping_OneMoreLocations.bam > $SAMPLELOC/$SAMPLE/${ACCESSION}_mapping_OneMoreLocations.sorted.bam
${SAM}samtools index $SAMPLELOC/$SAMPLE/${ACCESSION}_mapping_OneMoreLocations.sorted.bam

rm -rf $SAMPLELOC/$SAMPLE/Backmapping
rm -rf $SAMPLELOC/$SAMPLE/$SRA

# conduct PACVr Rscript
Rscript $PACVR -k $SAMPLELOC/$SAMPLE/${ACCESSION}.gb -b $SAMPLELOC/$SAMPLE/${ACCESSION}_mapping_OneMoreLocations.sorted.bam -t 0.5 -r TRUE -v TRUE -o $SAMPLELOC/$SAMPLE/${ACCESSION}_CoverageDepth.pdf

# merge metadata to one csv
${NCBIENT}esearch -db sra -query $SRA | efetch -format runinfo | cut -f20,7 -s -d, > $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/${ACCESSION}_SRA_meta.csv 
${NCBIENT}esearch -db nuccore -query $ACCESSION | efetch -format gb | sed -n '/Assembly Method/,/Sequencing/{/Sequencing/!p;}'  | awk '{$1=$1;print}' | sed 's/:: /\n/g' | paste -s -d ' ' | sed 's/  /\n/g' | awk /./ > $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/${ACCESSION}_assembly_tech.csv
find $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/ -type f -size 0 -exec rm {} \;
paste -d, $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/*.csv > $SAMPLELOC/$SAMPLE/${ACCESSION}.tmp/${ACCESSION}_metadata.csv
