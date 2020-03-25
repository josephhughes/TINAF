#!/bin/bash

directory=$1
ref=$2
bedpe=$3

for filepath in $directory/*R1_001.fastq.gz
do 
  echo $filepath
  file=$(basename $filepath)
  #echo $file
  stub=$(echo $file | cut -f1 -d "_")
  #echo $stub
  #echo TrimmomaticPE -phred33 $filepath ${filepath%R1_001.fastq.gz}R2_001.fastq.gz ${stub}_trimmomatic_R1.fastq ${stub}_unpaired1.fastq ${stub}_trimmomatic_R2.fastq ${stub}_unpaired2.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:25 MINLEN:5 -threads 8 
#  TrimmomaticPE -phred33 $filepath ${filepath%R1_001.fastq.gz}R2_001.fastq.gz ${stub}_trimmomatic_R1.fastq ${stub}_unpaired1.fastq ${stub}_trimmomatic_R2.fastq ${stub}_unpaired2.fastq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:25 MINLEN:5 -threads 8 

  #echo bwa mem $ref ${stub}_trimmomatic_R1.fastq ${stub}_trimmomatic_R2.fastq | samtools view -bS | samtools sort -o ${stub}_bwa.bam
#  bwa mem $ref ${stub}_trimmomatic_R1.fastq ${stub}_trimmomatic_R2.fastq | samtools view -bS | samtools sort -o ${stub}_bwa.bam
  
#  samtools index ${stub}_bwa.bam
  #echo /home1/orto01r/alicTest/bamclipper-master/bamclipper.sh -b ${stub}_bwa.bam -p $bedpe -s /software/samtools-1.9/bin/samtools
#  /home1/orto01r/alicTest/bamclipper-master/bamclipper.sh -b ${stub}_bwa.bam -p $bedpe -s /software/samtools-1.9/bin/samtools -n 8

  #echo weeSAM --bam ${stub}_bwa.bam --html ${stub}_bwa
#  weeSAM --bam ${stub}_bwa.primerclipped.bam --html ${stub}_bwa --overwrite
  
  #echo samtools view -F 0x04 -b ${stub}_bwa.bam > ${stub}_bwa_mapped.bam
#  samtools view -F 0x04 -b ${stub}_bwa.primerclipped.bam | samtools sort  -o  ${stub}_bwa_mapped.bam
  samtools index ${stub}_bwa_mapped.bam
  wysiwyg_consensus.pl -bam ${stub}_bwa_mapped.bam -ref $ref -stub ${stub}_cov5 -cov 5

done