#!/bin/bash
cat any_path/samples | while read sample
do
cutadapt \
--discard-untrimmed --no-indels  \
--overlap 15 \
-e 0.15 \
-g file:any_path/fw_primers \
-o any_path/demultiplexed/$sample.{name}.1.fastq.gz \
/Users/miguelcamacho/Library/CloudStorage/Dropbox/research/repos/tidyGenR/inst/extdata/raw/"$sample".1.fastq.gz
done > cutadapt.log
