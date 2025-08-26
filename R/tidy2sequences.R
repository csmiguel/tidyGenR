#' Output FASTA sequences
#'
#' Write FASTA and/or return BioString object with selected sequences.
#' @param td Tidy variants or genotypes. Dataframe with at least one column
#' with DNA sequences.
#' @param fasta_header Naming of FASTA headers. The string is passed to 'glue'
#'  for forming the FASTA headers and selecting distinct rows.
#' @param filename Path to write FASTA file, or FALSE, to prevent writing
#' DNA sequences to a FASTA file.
#' @param seq Name of the column in 'td' with DNA sequences.
#' @details
#' 'fasta_header' is a flexible selector and constructor for FASTA headers.
#' The variables included are used to filter distinct rows in 'td'.
#' Examples of outputs:
#'  - "\{locus\}_\{allele\}", all different allele sequences
#' per locus.
#'  - "\{sample\}_\{locus\}_\{variant\}", all sequences for all the
#'  variants for all samples, thus repeated sequences from the same variant
#'  corresponding to different samples.
#'  - "\{md5\}", all different DNA sequences.
#'  - "\{sample\}", one sequence per sample. Since one sample matches
#'  to many sequences, the first occurrence in the dataframe is selected.
#' @return 'DNAStringSet' object with selected  sequences.
#' @examples
#' data("genotypes")
#' tidy2sequences(
#'     td = genotypes,
#'     fasta_header = "{locus}_{allele}",
#'     filename = FALSE
#' )
#' @export
tidy2sequences <- function(td,
                           fasta_header = "{locus}_{allele}",
                           filename = FALSE,
                           seq = "sequence") {
    stopifnot(
        "data.frame" %in% class(td),
        seq %in% names(td)
    )
    td <- drop_na(td)
    # filter tidy data according to output desired
    # variables to select
    col_select <-
        str_remove_all(
            str_extract_all(
                string = fasta_header,
                pattern = "\\{[a-zA-Z0-9]*\\}"
            )[[1]],
            "\\{|\\}"
        )
    # filter distinct elements
    td_distinct <-
        distinct(td,
            across(all_of(col_select)),
            .keep_all = TRUE
        )
    # get names for fasta headers
    names_fh <-
        as.character(
            mutate(td_distinct,
                fh = glue(fasta_header)
            )[["fh"]]
        )
    # create DNAStringSet
    bs <-
        DNAStringSet(pull(td_distinct, seq))
    names(bs) <- names_fh
    if (is.character(filename)) {
        # write fasta
        writeXStringSet(bs,
            format = "fasta",
            filepath = filename
        )
        message("Sequences have been written to: ", filename)
    } else if (!filename) {
        message("Sequences have not been written to fasta file.")
    }
    return(bs)
}
