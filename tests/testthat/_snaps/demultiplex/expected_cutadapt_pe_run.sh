#!/bin/bash
cat any_path/samples | while read sample
do
cutadapt \
--discard-untrimmed --pair-adapters --no-indels  \
--overlap 15 \
-e 0.15 \
-g file:any_path/fw_primers \
-G file:any_path/rv_primers \
-o any_path/demultiplexed_pe_run/$sample.{name}.1.fastq.gz \
-p any_path/demultiplexed_pe_run/$sample.{name}.2.fastq.gz \
/Users/miguelcamacho/Library/CloudStorage/Dropbox/research/repos/tidyGenR/inst/extdata/raw/"$sample".1.fastq.gz \
/Users/miguelcamacho/Library/CloudStorage/Dropbox/research/repos/tidyGenR/inst/extdata/raw/"$sample".2.fastq.gz
done > any_path/cutadapt_pe_run.log
