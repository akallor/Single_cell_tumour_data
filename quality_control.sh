#!/usr/bin/env bash
#Perform quality control and read trimming on the reads
#Using trim-galore but also try with FASTQC

trim_galore --paired --cores 4 \
  Tumor1_R1.fastq.gz Tumor1_R2.fastq.gz \
  -o trimmed/
