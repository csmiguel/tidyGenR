#' Primers for multiplex PCR
#'
#' A dataframe containing 27 primer pairs used to amplify intron loci in
#' *Rattus baluensis* (Igea et al. 2010; Camacho-Sanchez et al. 2018).
#' @references
#' Camacho-Sanchez et al. (2018). _Interglacial refugia on tropical
#' mountains: novel insights from the summit rat (<i>Rattus baluensis</i>),
#'  a Borneo mountain endemic_. Diversity and Distributions.
#'
#' Igea et al. (2010) _Novel intron markers to study the phylogeny of
#' closely related mammalian species_. BMC Evolutionary Biology.
#' @format A data frame with 27 rows and 3 variables.
#' \describe{
#'   \item{locus}{Locus name.}
#'   \item{fw}{5' to 3' forward primer sequence.}
#'   \item{rv}{5' to 3' reverse primer sequence.}
#' }
"primers"

#' Variants
#'
#' Filtered tidy variants of *Rattus baluensis* from Trusmadi
#' (Camacho-Sanchez et al. _in preparation_),
#' corresponding to 27 loci from 44 samples.
#' @format tidy dataframe or tibble.
#' \describe{
#'   \item{sample}{Sample name.}
#'   \item{locus}{Locus name.}
#'   \item{variant}{Variant name.}
#'   \item{reads}{Number of reads supporting allele.}
#'   \item{nt}{Length in number of nucleotides of allele sequence.}
#'   \item{md5}{MD5 checksum of allele sequence.}
#'   \item{sequence}{DNA sequence of allele.}
#' }
"variants"

#' Per-locus truncation lengths for forward and reverse reads.
#'
#' Truncation lengths for forward and reverse reads for
#' the dataset containing 30 primer pairs used to amplify intron loci in
#' *Rattus baluensis* (Camacho-Sanchez et al. _in preparation_).
#'
#' @format Data frame with three columns:
#' \describe{
#'   \item{locus}{Locus name.}
#'   \item{trunc_f}{Truncation length for forward reads.}
#'   \item{trunc_r}{Truncation length for reverse reads.}
#' }
"trunc_fr"

#' Genotypes
#'
#' Filtered genotypes of 20 loci and 91 samples of *Rattus baluensis* used in
#'  Camacho-Sanchez et al. _in preparation_.
#'
#' @format tidy dataframe or tibble.
#' \describe{
#'   \item{sample}{Sample name.}
#'   \item{locus}{Locus name.}
#'   \item{allele}{Allele name.}
#'   \item{allele_no}{Allele ordinal number.}
#'   \item{reads}{Number of reads supporting allele.}
#'   \item{nt}{Length in number of nucleotides of allele sequence.}
#'   \item{md5}{MD5 hash of allele sequence.}
#'   \item{sequence}{DNA sequence of allele.}
#' }
"genotypes"

#' List of dataframes with _qvalues_.
#'
#' List of 6 dataframes with 91 samples containing q-ancestry values
#' from STRUCTURE runs.
#' @format A list of 6 dataframes with 91 rows.
"qlist"

#' List of tidy variants
#'
#' List of 5 tidy variants dataframes with variant calls using
#' different strategies.
#' @format A list of 5 dataframes.
#'  \describe{
#'   \item{ind_bs16}{pool = F; band_size = 16}
#'   \item{ind_bs0}{pool = F; band_size = 0}
#'   \item{pool_bs16}{pool = T; band_size = 16}
#'   \item{pool_bs0}{pool = T; band_size = 0}
#'   \item{amplisas}{AmpliSAT}
#' }
"variant_calls"

#' MSA from locus tmem87a
#'
#' Multiple sequence alignment stored as Biostrings 'DNAStringSet'
#' from locus tmem87a from _Rattus baluensis_.
#' @format DNAStringSet with 179 aligned sequences.
"msa"
