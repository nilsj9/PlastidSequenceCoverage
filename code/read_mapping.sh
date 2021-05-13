#!/bin/bash

# REQUIREMENTS
#ncbi-entrez-direct 
#sra-toolkit >2.10.8 
#samtools >1.10
#BOWTIE2 >2.3.4.1
#TRIMMOMATIC

# SET PATHS
SAM=/PATH/TO/SAMTOOLS
SRAT=/PATH/TO/SRA-TOOLKIT
NCBIENT=/PATH/TO/NCBI-ENTREZ-DIRECT
BOWTIE=/PATH/TO/BOWTIE2
TRIMMOMATIC=/PATH/TO/TRIMMOMATIC
SAMPLELOC=/PATH/TO/SAMPLES

# DEFINITIONS
SAMPLE=$1
ACCESSION=$2
SRA=$3

# FOLDERS
mkdir -p $SAMPLELOC/$SAMPLE/ReadMapping

# download GenBank
${NCBIENT}esearch -db nucleotide -query $ACCESSION | ${NCBIENT}efetch -format gb > $SAMPLELOC/$SAMPLE/${ACCESSION}.gb

# download reference fasta
${NCBIENT}esearch -db nucleotide -query $ACCESSION | ${NCBIENT}efetch -format fasta > $SAMPLELOC/$SAMPLE/${ACCESSION}_ref.fasta

# download reads
cd $SAMPLELOC/$SAMPLE
${SRAT}prefetch --max-size 50000000 $SRA
cd /scratch/nilsj24/
${SRAT}fasterq-dump.2.10.8 --split-3 --skip-technical $SAMPLELOC/$SAMPLE/$SRA/$SRA.sra -O $SAMPLELOC/$SAMPLE/ReadMapping

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}.sra_1.fastq $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}.sra_2.fastq  $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_1_PE.fastq $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_1_SE.fastq $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_2_PE.fastq $SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_2_SE.fastq ILLUMINACLIP:/$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

# conduct mapping of reads
REF=$SAMPLELOC/$SAMPLE/${ACCESSION}_ref.fasta
READSPER1=$SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_1_PE.fastq
READSPER2=$SAMPLELOC/$SAMPLE/ReadMapping/${SRA}_2_PE.fastq

mkdir -p $SAMPLELOC/$SAMPLE/ReadMapping/db
${BOWTIE}bowtie2-build $REF $SAMPLELOC/$SAMPLE/ReadMapping/db/$ACCESSION
${BOWTIE}bowtie2 -x $SAMPLELOC/$SAMPLE/ReadMapping/db/$ACCESSION -1  $READSPER1 -2 $READSPER2 -S $SAMPLELOC/$SAMPLE/ReadMapping/${ACCESSION}_mapping.sam
${SAM}samtools view -Sb -F 0x04 $SAMPLELOC/$SAMPLE/ReadMapping/${ACCESSION}_mapping.sam > $SAMPLELOC/$SAMPLE/ReadMapping/${ACCESSION}_mapping_OneMoreLocations.bam
${SAM}samtools sort $SAMPLELOC/$SAMPLE/ReadMapping/${ACCESSION}_mapping_OneMoreLocations.bam > $SAMPLELOC/$SAMPLE/${ACCESSION}_mapping_OneMoreLocations.sorted.bam
${SAM}samtools index $SAMPLELOC/$SAMPLE/${ACCESSION}_mapping_OneMoreLocations.sorted.bam

rm -rf $SAMPLELOC/$SAMPLE/ReadMapping
rm -rf $SAMPLELOC/$SAMPLE/$SRA
