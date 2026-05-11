#!/usr/bin/env bash
#Download SRA Toolkit (ignore if already downloaded)
#After downloading SRA toolkit:
prefetch SRS18169578 SRS18169579 SRS18169580
fasterq-dump --split-files SRS18169578
fasterq-dump --split-files SRS18169579
fasterq-dump --split-files SRS18169580
pigz *.fastq   # compress in parallel
#TODO: Download the Sus scrofa reference: Ensembl Sscrofa11.1 (GCA_000003025.6), including the GTF annotation file.
