*Pipeline for Inferring Plastid Sequencing Depth and Evenness*
==============================================================

A bioinformatic pipeline for inferring sequencing depth and evenness on complete plastid genomes


Description of Pipeline
-----------------------

This pipeline illustrates the individual steps taken during the analysis of sequencing depth and evenness of a set of 194 genome records of complete plastid genomes. The pipeline comprises of four parts: (i) read mapping, (ii) metadata extraction, (iii) calculation of sequencing depth and evenness, and (iv) statistical analysis.


Part i: Read mapping
--------------------

Input: _samples\_list.csv_

##### Conducting read mapping on a series of plastid genome records
```
#!/bin/bash

INPUT=/PATH/TO/samples_list.csv
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read -r SRA SAMPLE ACCESSION;
do
bash ./read_mapping.sh $SAMPLE $ACCESSION $SRA </dev/null
done < $INPUT
IFS=$OLDIFS
```

Part ii: Metadata extraction
----------------------------

Input: _samples\_list.csv_

##### Extract metadata of a series of plastid genome records
```
#!/bin/bash

INPUT=/PATH/TO/samples_list.csv
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read -r SRA SAMPLE ACCESSION;
do
bash ./metadata_extraction.sh $SAMPLE $ACCESSION $SRA </dev/null
done < $INPUT
IFS=$OLDIFS
```


Part iii: Sequencing depth and evenness
---------------------------------------

Input: _Customize path to resources in code_

##### Gather all information into a single csv file.
```
Rscript depth_and_evenness.R
```


Part iv: Statistical analysis
-----------------------------

Input: _stats\_ready\_metadata.csv_

##### Conduct statistical analysis on a series of plastid genome records
```
Rscript plastid_analysis.R
```
