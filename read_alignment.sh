#!/usr/bin/env bash

#Perform read alignment post-quality trimming and QC using STAR

STAR --runMode genomeGenerate \
  --genomeDir sscrofa_index/ \
  --genomeFastaFiles Sscrofa11.1.dna.fa \
  --sjdbGTFfile Sus_scrofa.Sscrofa11.1.gtf

# STARsolo alignment (10x v3 example)
STAR --soloType CB_UMI_Simple \
  --soloCBwhitelist 3M-february-2018.txt \
  --readFilesIn Tumor1_R2.fastq.gz Tumor1_R1.fastq.gz \
  --genomeDir sscrofa_index/ \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix Tumor1/
