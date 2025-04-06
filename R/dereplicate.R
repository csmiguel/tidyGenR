#' De-replicate reads into frequency tables for a set of FASTQ files
#'
#' Reads FASTQ files and computes frequency of unique sequences per sample.
#' Depth for unique sequences are organized in *samples* x *unique sequences* tables.
#' Unique Ssequences are ordered in descending frequency.
#' @param fs Character vector with paths to all FASTQ to de-replicate.
#' @param min_sam_fr Numeric. Minimum number of sequence counts in a sample to
#' be retained (cell number).
#' @param min_loc_fr Numeric. Minimum frequency of de-replicated sequence in the group to be retained.
#' If \eqn{min\_loc\_fr \in (0,1)}, then a proportion relative to the
#' most frequent sequence is applied.
#' @param by Regex pattern to group FASTQ files in the list. Passed to `stringr::str_extract()`.
#' @param out_xlsx File name to write tables with de-replicated sequences
#' (*Default: NULL; no file is written*).
#' @details
#' The *by* parameter allows flexible grouping of files in the list. However, the results are not added within each group; individual results for each sample are always returned.
#' For instance, given 3 files s1_loc1_F.fq, s1_loc1_R.fq and s2_loc1_F.fq:
#' - \code{"([a-zA-Z0-9]*_[a-zA-Z0-9]*)"}, returns *s1_loc1* and *s2_loc2*.
#' - \code{"_([a-zA-Z0-9]*)_"}, returns *loc1* and *loc2*.
#' - \code{"([a-zA-Z0-9]*_[F|R])"}, returns *loc1_F*, *loc1_R* and *loc2_F*.
#' The `min_sam_fr` and `min_loc_fr` filters drop data not passing the filters, so they will become
#' zero or absent when combined.
#' If a path to an EXCEL file is set, each element in the list is written to a
#' different sheet in the workbook.
#' @returns List of extracted groups (see 'by'). Each element is the list is
#' a dataframe with:
#'  - column 1: 'sequence', DNA sequence of read.
#'  - column 2: 'md5', md5 hash of DNA sequence.
#'  - column >= 3: frequency (*integers*) of sequences per sample passing
#' 'min_sam_fr' and 'min_loc_fr' filters.
#' @examples
#' fq <-
#'  list.files(system.file("extdata", "truncated",
#'                         package = "tidyGenR"),
#'                         pattern = "fastq.gz",
#'             full.names = TRUE)
#' dereplicate(fq)
#' @export
dereplicate <- function(
    fs, min_sam_fr = 2,
    min_loc_fr = 0.001,
    by = "_([a-zA-Z0-9]*_[F|R])",
    out_xlsx = NULL) {
    # get groups
    groups_dereps <-
        unique(str_extract(
            basename(fs),
            by, 1
        ))
    # dataframe for each group with dereplicated reads
    list_dereps <-
        lapply(groups_dereps, function(y) {
            fsel <-
                grep(y, fs, value = TRUE)
            dr <- derepFastq(fsel) |>
                dada2list(basename(fsel))
            un <-
                lapply(names(dr), function(x) {
                    z <- dr[[x]]$uniques
                    data.frame(
                        sample = x,
                        sequence = names(z),
                        freq = z
                    )
                })
            w <-
                do.call(what = "rbind", un) |>
                filter(.data$freq >= min_sam_fr) |>
                pivot_wider(
                    id_cols = "sequence",
                    names_from = "sample",
                    values_from = "freq",
                    values_fill = 0
                )
            # reorder sequences according to the sum of counts
            ww <- w[order(rowSums(w[, -1]), decreasing = TRUE), ]
            rs <- rowSums(ww[, -1])
            if (min_loc_fr < 1 & min_loc_fr > 0) {
                ww <- ww[rs / max(rs) >= min_loc_fr, ]
            } else if (min_loc_fr > 1) {
                ww <- ww[rs >= min_loc_fr, ]
            }
            ww <-
                rowwise(ww) |>
                mutate(md5 = digest(.data$sequence), .after = .data$sequence)
            return(ww)
        })
    names(list_dereps) <- groups_dereps
    if (!is.null(out_xlsx)) {
        write_xlsx(list_dereps, path = out_xlsx)
        message(
            "De-replicated sequences written to ", out_xlsx, ".\n",
            "If the number of sheets in too large EXCEL might crash at",
            " opening."
        )
    }
    return(list_dereps)
}
