#!/bin/bash
for bamfile in *simualted_hg19.clean.bam
do
  
  NAME=${bamfile}
  echo $NAME  
# get just R1:
  samtools view  -@ 1 -f 0x40 -O "BAM" -o ${NAME}_R1.bam $NAME
  # turn this to bed:
  bamToBed -ed -i ${NAME}_R1.bam > ${NAME}_R1.bed

done

