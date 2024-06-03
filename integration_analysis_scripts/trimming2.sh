#!/bin/bash

sample1=$1
sample2=$2
NAME=$3
echo $sample1
echo $sample2
# do the trimming and control operations

#fixing read names
/common/help_software/BCBIO/tools/bin/reformat.sh -Xmx60g in=$sample1 in2=$sample2 out=${NAME}_1_fixname.fastq.gz out2=${NAME}_2_fixname.fastq.gz trimreaddescription=t addslash=t slashspace=f
#for reads 1, remove everything left of viral sequences and everything right of linker sequence, then remove ends with quality < 10:
/common/help_software/BCBIO/tools/bin/bbduk.sh -Xmx60g in=${NAME}_1_fixname.fastq.gz out=${NAME}_1_decontam1.fastq.gz ref=NC_001802.1.fasta ktrim=l 
/common/help_software/BCBIO/tools/bin/bbduk.sh -Xmx60g in=${NAME}_1_decontam1.fastq.gz out=${NAME}_1_decontam2.fastq.gz ref=linker.fasta ktrim=r
/common/help_software/BCBIO/tools/bin/bbduk.sh -Xmx60g in=${NAME}_1_decontam2.fastq.gz out=${NAME}_1_decontam_trimmed.fastq.gz qtrim=rl trimq=10 threads=10
# for reads 2, remove everything left of linker and everything right of viral sequence, then remove ends with quality < 10:
/common/help_software/BCBIO/tools/bin/bbduk.sh -Xmx60g in=${NAME}_2_fixname.fastq.gz out=${NAME}_2_decontam1.fastq.gz ref=NC_001802.1.fasta ktrim=r mink=9
/common/help_software/BCBIO/tools/bin/bbduk.sh -Xmx60g in=${NAME}_2_decontam1.fastq.gz out=${NAME}_2_decontam2.fastq.gz ref=linker.fasta ktrim=l mink=9
/common/help_software/BCBIO/tools/bin/bbduk.sh -Xmx60g in=${NAME}_2_decontam2.fastq.gz out=${NAME}_2_decontam_trimmed.fastq.gz qtrim=rl trimq=10 threads=10
# re-pair reads
/common/help_software/BCBIO/tools/bin/repair.sh -Xmx60g in=${NAME}_1_decontam_trimmed.fastq.gz in2=${NAME}_2_decontam_trimmed.fastq.gz out=${NAME}_1_clean.fastq.gz out2=${NAME}_2_clean.fastq.gz repair=t
