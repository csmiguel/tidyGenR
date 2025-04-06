#' Convert sequence-based variants to length-based
#'
#' @details Recodes dataframe with sequence-based variants, such as the output
#' from 'variant_call' to
#' length-based variants. The sequence length 'nt', instead of the nucleotide
#'  sequence, is used for determining variants.
#' Each variant is named after its number of nucleotides.
#'  Number of reads from previous variants with equal lengths are aggreated.
#' Legacy for microsatellite data.
#' @param variants Dataframe with tidy variants.
#' @examples
#' data("variants")
#' var_seq2len(variants)
#' @export
#' @return Dataframe with length-based variants and the variables:
#' 'locus' 'sample' 'variant' 'reads'.
var_seq2len <- function(variants) {
    mand_vars <- c("locus", "sample", "reads", "nt")
    if (!all(mand_vars %in% names(variants))) {
        stop(
            "'variant' names must include '",
            paste(mand_vars, collapse = " "), "'"
        )
    }
    a <-
        select(variants, .data$locus, .data$sample, .data$reads, .data$nt) |>
        group_by(.data$locus, .data$sample, .data$nt) |>
        summarize(reads = sum(.data$reads)) |>
        rename(variant = .data$nt) |>
        ungroup()
    return(a)
}
