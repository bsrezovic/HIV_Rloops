#!/bin/bash
for bamfile in *_hg19.clean.bam
do
  
  NAME=${bamfile%.bam}
  echo $NAME  
# get just R1:
  samtools view  -@ 1 -f 0x40 -O "BAM" -o ${NAME}_R1.bam ${NAME}.bam
  # turn this to bed:
  bamToBed -ed -i ${NAME}_R1.bam > ${NAME}_R1.bed

done

