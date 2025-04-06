#' Output reads across pipeline
#'
#' It generates a dataframe with tracked reads per samples per locus.
#' through the pipeline, from dataframes or FASTQ files.
#'
#' @param l List of elements to retrieved reads from.
#' @param sample_pattern Pattern used to extract sample names from FASTQ files.
#' @return Dataframe with tracked reads through the pipeline. The first column
#' is named 'sample' and the following
#' @details 'l' is a named list. It can be either a dataframe with two
#' mandatory fields ('sample' and 'reads'),
#' or a character vector with 2 elements: 1, the path to a directory containing
#' FASTQ files; and
#' ,2, a pattern matching desired FASTQ files. By, default, sample names
#'  matching
#' 'sample_pattern'  '^\[A-Za-z0-9\]*' are extracted from FASTQ files.
#' 'left_join()' is used
#' iteratively from the first to the last element in 'l', so the elements in 'l'
#' need to be nested from element 1 to the last.
#' @examples
#' # from folder with FASTQ files
#' path2trunated <-
#'    system.file("extdata", "truncated",
#'              package = "tidyGenR")
#' l <- list(
#'       truncated = c(path2trunated, "F_filt.fastq.gz"))
#' reads_track(l)
#' @export
reads_track <- function(l, sample_pattern = "^[A-Za-z0-9]*") {
    l2 <-
        lapply(l, function(x) {
            if (inherits(x, "character")) {
                fp <- list.files(x[1], x[2], full.names = TRUE)
                reads <-
                    countFastq(fp)[, "records", drop = FALSE] |>
                    rownames_to_column("sample") |>
                    mutate(sample = str_extract(
                        sample,
                        sample_pattern
                    )) |>
                    group_by(sample) |>
                    summarise(temp = sum(.data$records))
                names(reads)[2] <- x[1]
            } else if (inherits(x, "data.frame")) {
                reads <-
                    group_by(x, sample) |>
                    summarise(variants = sum(reads))
            }
            return(reads)
        })

    z <- l2[[1]]
    if (length(l) > 1) {
        for (i in seq_len(length(l2) - 1)) {
            z <- left_join(z, l2[[i + 1]], by = "sample")
        }
    }
    names(z) <- c("sample", names(l))
    return(z)
}
