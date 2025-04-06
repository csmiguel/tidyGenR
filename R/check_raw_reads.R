#' Check input raw FASTA/Q files
#'
#' @details
#' Multiple checks on raw input FASTQ:
#'  - Reads in FASTQ above threshold.
#'  - Unique sample names.
#'  - Not orphan F/R files for paired reads.
#'  For single-end data check 3 is excluded and NULL is returned in `checks[4]`.
#' @param freads Character vector with file paths to forward reads.
#' @param rreads Character vector with file paths to reverse reads.
#' @param low_readcount Numeric threshold to warn on number of reads in FASTA
#'  (DEFAULT = 10).
#' @return List with 4 elements:
#'  - 'f_reads', character vector with basename of forward reads.
#'  - 'r_reads', character vector with basename of reads reads
#'  or NULL for single-end experiments.
#'  - 'snames', character vector with sample names.
#'  - 'checks', Logical vector with passed checks 1-3.
#' @examples
#' freads <-
#' list.files(system.file("extdata", "raw", package = "tidyGenR"),
#' pattern = "1.fastq.gz",
#'     full.names = TRUE)
#' check_raw_reads(freads)
#' @export
check_raw_reads <- function(freads, rreads = NULL, low_readcount = 10) {
    # list with file names.
    a <- list(freads = freads, rreads = rreads)
    stopifnot(!all(vapply(a, is.null, logical(1))))
    # list with basenames.
    reads_bs <-
        lapply(a, function(x) {
            if (!is.null(x)) {
                basename(x)
            } else {
                NULL
            }
        })
    # list with sample names
    snames <-
        lapply(reads_bs, function(x) {
            if (!is.null(x)) {
                gsub("(^[a-zA-Z0-9]*).*$", "\\1", x)
            } else {
                NULL
            }
        })
    # check 1: number of reads
    ch1 <- check1_lowcount(a, low_readcount)
    # check 2: duplicated sample names
    ch2 <- check2_unique_snames(snames)
    # check 3: all F files have their corresponding R file.
    ch3 <- check3_non_orphan(a, snames)
    # return files and sample names
    slist <-
        list(
            f_reads = reads_bs$freads,
            r_reads = reads_bs$rreads,
            samples = snames[[1]],
            checks = c(check1 = ch1, check2 = ch2, check3 = ch3)
        )
    return(slist)
}


#' Check1
#'
#' All FASTQ have a number of reads above threshold.
#' @param list_fp List of length 1 or 2 character vectors with
#' absolute F (and R) file paths.
#' @param low_readcount 'low_readcount' from 'check_raw_reads()'.
#' @return Logical, TRUE if test is passed. Warning if FALSE.
check1_lowcount <- function(list_fp, low_readcount = low_readcount) {
    probe <-
        lapply(list_fp, function(x) {
            if (!is.null(x)) {
                all(countFastq(x)$record > low_readcount)
            } else {
                NULL
            }
        })
    if (all(unlist(probe))) {
        message(
            "All F and (R files) passed check on number of reads above ",
            low_readcount, "."
        )
        return(TRUE)
    } else {
        warning(
            "Some files have less than ",
            low_readcount, " reads."
        )
        return(FALSE)
    }
}

#' Check2
#'
#' Sample names are unique.
#' @param snames List of length 1 or 2 character vectors with sample names.
#' @return Logical, TRUE if test is passed. STOPs if FALSE.
check2_unique_snames <- function(snames = snames) {
    probe <-
        vapply(
            snames, function(x) {
                any(duplicated(x))
            },
            logical(1)
        )
    if (all(!probe)) {
        message("Sample names are unique.")
        return(TRUE)
    } else {
        stop(
            "Duplicated sample names:",
            snames[[1]][duplicated(snames[[1]])],
            snames[[2]][duplicated(snames[[2]])]
        )
    }
}

#' Check3
#'
#' Not orphan F/R FASTQ files.
#' @param list_fp List of length 1 or 2 character vectors with
#' absolute F (and R) file paths.
#' @param snames List of length 1 or 2 character vectors with sample names.
#' @return Logical, TRUE if test is passed. Warning if FALSE.
#' NULL for single end reads.
check3_non_orphan <- function(list_fp, snames = snames) {
    probe <- any(vapply(list_fp, is.null, logical(1)))
    if (isTRUE(probe)) {
        message("Only F reads in check 3: orphan reads cannot be evaluated.")
        return(NULL)
    } else {
        probe2 <- identical(sort(snames[[1]]), sort(snames[[2]]))
        if (probe2) {
            message("All F files have their corresponding R file.")
            return(TRUE)
        } else {
            warning("NOT all F files have their corresponding R file.")
            return(FALSE)
        }
    }
}
