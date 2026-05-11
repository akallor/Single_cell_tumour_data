#!/usr/bin/env bash
#Download SRA Toolkit (ignore if already downloaded)
#After downloading SRA toolkit:
prefetch SRS18169578 SRS18169579 SRS18169580
fasterq-dump --split-files SRS18169578
fasterq-dump --split-files SRS18169579
fasterq-dump --split-files SRS18169580
pigz *.fastq   # compress in parallel
