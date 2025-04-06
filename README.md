[![R-CMD-check](https://github.com/csmiguel/tidyGenR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/csmiguel/tidyGenR/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/csmiguel/tidyGenR/branch/main/graph/badge.svg)](https://codecov.io/gh/csmiguel/tidyGenR)
[![License: GPL v3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<img src="man/figures/tidyGenR.png" alt="hex sticker" style="width:20%;"/>

# tidyGenR: Tidy multilocus amplicon genotypes in R

R package optimized for genotyping multilocus amplicon sequences from High Throughput Sequencing in diploid organism. Starting with FASTQ multilocus amplicon sequencing reads from high throughput sequencing, it demultiplexes reads by locus, call variants (using DADA2), and genotypes samples. It generates tidy variants and genotypes outputs, and it provides a family of functions for data manipulation and format conversion.

* optimized to genotype amplicon reads from diploid markers.
* FASTQ inputs: overlapping and non-overlapping, single-end or paired-end reads from Illumina.
* genotyping of one or multiple species simoultaneously
* various output format: tidy genotypes, FASTA, STRUCTURE, Genalex.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("csmiguel/tidyGenR")
```

# Workflow

<img src="man/figures/fig1_flow.svg" alt="workflow diagram" style="width:100%;"/>

