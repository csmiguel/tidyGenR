#' Align variant sequences against a reference FASTA database
#'
#' Performs pairwise alignments between DNA sequence variants contained in a
#' tidyGenR variants table and a reference FASTA database using
#' \code{DECIPHER::AlignPairs()}. The best matching reference sequence(s) for
#' each variant are returned and appended as additional columns to the variants
#' table.
#'
#' This function is intended for exploratory annotation of sequence variants,
#' such as matching alleles against reference haplotypes, species databases, or
#' known locus variants.
#'
#' @param variants A tidyGenR variants table containing DNA sequence variants.
#' @param ref_fasta Path to a FASTA file containing reference DNA sequences.
#' @param n_best Integer specifying the number of top-scoring alignments to
#'   retain for each query sequence. Defaults to \code{1}.
#' @param reduced Logical. If \code{TRUE} (default), only a reduced set of
#'   alignment statistics is returned (\code{query}, \code{ref},
#'   \code{Score}, \code{Matches}, and \code{Mismatches}).
#'   If \code{FALSE}, all columns returned by
#'   \code{DECIPHER::AlignPairs()} are retained.
#' @param ... Additional arguments passed to
#'   \code{DECIPHER::AlignPairs()}.
#'
#' @details
#' Reference sequences are read using
#' \code{Biostrings::readDNAStringSet()}. Variant sequences are extracted from
#' the tidyGenR variants object using \code{tidy2sequences()} and aligned
#' against all reference sequences using pairwise alignment.
#'
#' Query sequences are internally identified using MD5 hashes to facilitate
#' unambiguous tracking across the workflow.
#'
#' @return
#' A tidy variants table with additional columns containing alignment
#' information for the best matching reference sequence(s).
#'
#' @examples
#'data("variants")
#'ref_al <-
#'     system.file("extdata/reference_alleles.fasta", package = "tidyGenR")
#'align_variants_ref(variants, ref_al)
#'
#' @export
align_variants_ref <- function(variants, ref_fasta, n_best = 1, reduced = TRUE, ...) {
  # read fasta ref seqs
  ref_seqs <-
    readDNAStringSet(ref_fasta)
  stopifnot(length(ref_seqs) > 0)
  # control statements
  if(length(ref_seqs) < n_best)
    stop("The number of ref sequences needs to be larger than 'n_best'")
  if(n_best > 5)
    warning("A maximum of n_best = 3 is recommended")
  # get query sequences from variants object
  query_seqs <-
    suppressMessages(tidy2sequences(variants, fasta_header = "{md5}"))
  stopifnot(length(query_seqs) > 0)
  # index for pw alignment
  df_index <-
    cross_join(
      data.frame(Pattern = seq_along(ref_seqs)), # reference
      data.frame(Subject = seq_along(query_seqs)), # query
    )
  # run alignment
  al <-
    AlignPairs(ref_seqs,
               query_seqs,
               pairs = df_index,
               ...)

  # rename Pattern (ref) and Subject (query) sequences from index [0-9]+ to their sequence names
  al_renamed <-
    al |>
    mutate(
      ref = mapvalues(Pattern,
                      from = seq_along(ref_seqs),
                      names(ref_seqs)),
      query = mapvalues(Subject,
                        from = seq_along(query_seqs),
                        names(query_seqs)))

  # get the best alignment for each query sequence
  al_best <-
    al_renamed |>
    ddply(~Subject, function(x) {
      x[order(x$Score, decreasing = TRUE),][seq_len(n_best),]
    })
  if(reduced)
    al_best <-
    select(al_best, query, ref, Score, Matches, Mismatches)

  # join to variants
  vars_al <-
    variants |>
    left_join(al_best,
              by = c("md5" = "query"),
              relationship = "many-to-many")


  return(vars_al)
}
