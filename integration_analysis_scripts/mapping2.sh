#!/bin/bash

#sample1=$1
#sample2=$2
NAME=$1
GENOME=/common/genomes/data/hg19/hg19.fa
HIV=NC_001802.1.fasta
PICARD=/common/Alignment_software/picard/build/libs/picard.jar
#unipping and changing the names because bwa requires unzipped files with .fq extension:
gunzip -k ${NAME}_1_clean.fastq.gz
mv ${NAME}_1_clean.fastq ${NAME}_R1.fq
gunzip -k ${NAME}_2_clean.fastq.gz
mv ${NAME}_2_clean.fastq ${NAME}_R2.fq

#Map to HG19
bwa mem -t 22 -M $GENOME ${NAME}_R1.fq ${NAME}_R2.fq -o ${NAME}.sam
samtools view ${NAME}.sam -u -@ 10 -o ${NAME}_2.sam 
samtools sort ${NAME}_2.sam -O BAM -m 5G -@ 10 -o ${NAME}_hg19.sorted.bam
/common/help_software/sambamba0.6.1/sambamba_v0.6.1 view -F "mapping_quality >= 30 and proper_pair and template_length < 900 and not (unmapped or secondary_alignment or supplementary)" -f bam -t 10 ${NAME}_hg19.sorted.bam -o ${NAME}_hg19.filtered.bam

#map to HIV as well
bwa mem -t 22 -M $HIV ${NAME}_R1.fq ${NAME}_R2.fq -o ${NAME}_HIV1.sam 
samtools view ${NAME}_HIV1.sam -u -@ 10 -o ${NAME}_HIV1_2.sam 
samtools sort ${NAME}_HIV1_2.sam -O BAM -m 5G -@ 10 -o ${NAME}_HIV1.sorted.bam
/common/help_software/sambamba0.6.1/sambamba_v0.6.1 view -F "mapping_quality >= 30 and not (unmapped or secondary_alignment or supplementary)" -f bam -t 10 ${NAME}_HIV1.sorted.bam -o ${NAME}_HIV1.filtered.bam
samtools view ${NAME}_HIV1.filtered.bam | awk '{if(length($10)>20)print $1}' | sort | uniq > ${NAME}_HIV1_reads.txt

# remove reads obtained in previous step from hg19 mapping:
NRREADS=`wc -l < ${NAME}_HIV1_reads.txt`
echo "Number of reads is"
echo $NRREADS
if [[ $NRREADS -gt 0 ]]
    then
      java -jar $PICARD FilterSamReads -INPUT ${NAME}_hg19.filtered.bam -OUTPUT ${NAME}_hg19.nohivreads.bam -READ_LIST_FILE ${NAME}_HIV1_reads.txt -FILTER excludeReadList -TMP_DIR ./tmp
    else
      mv ${NAME}_hg19.filtered.bam ${NAME}_hg19.nohivreads.bam
fi
  
# deduplicate:
/common/help_software/sambamba0.6.1/sambamba_v0.6.1 markdup -r -t 10 --tmpdir=tmp/ ${NAME}_hg19.nohivreads.bam ${NAME}_hg19.clean.bam
