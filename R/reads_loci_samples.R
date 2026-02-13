#' Output reads per locus
#'
#' Creates dataframe with reads distribution across samples (rows) and loci (columns).
#' @param path Folder with FASTQ files or tidy variants.
#' @param pattern_fq Pattern to select FASTQ files from the folder
#' @param sample_locus Sample and locus patterns to extract from file names.
#' Regex group 1 is assumed to be sample name and group 2 locus name.
#' By default, both groups are a alphanumeric string of any length separated by
#' any non-alphanumeric character.
#' @param all.variants (T/F) TRUE, one column per variant;
#' F, one column per locus. Only for tidy variants.
#' @param var_id 'c('allele', 'md5', 'sequence')' column to use for naming
#' variant.
#' @return Dataframe with columns being loci and rows being samples.
#' Samples are in first column.
#' @examples
#' reads_loci_samples(
#'      path = system.file("extdata", "truncated",
#'                        package = "tidyGenR"),
#'      pattern_fq = "F_filt.fastq.gz")
#' # with variants dataframe
#' data("variants")
#' reads_loci_samples(path = variants)
#'
#' @export
reads_loci_samples <- function(path,
                               pattern_fq = "1.fastq.gz",
                               sample_locus = "(^[a-zA-Z0-9]*).([a-zA-Z0-9]*)",
                               all.variants = FALSE,
                               var_id = "md5") {
    if (inherits(path, "character")) {
        fp <- list.files(path, pattern_fq, full.names = TRUE)
        if (length(fp) == 0) {
            stop("No files found.")
        }
        message(
            "'all.variants' and 'var_id' are ignored when setting ",
            "'path' to a folder."
        )
        reads <-
            countFastq(fp)[, "records", drop = FALSE] |>
            rownames_to_column("a") |>
            mutate(
                sample = str_extract(.data$a,
                    sample_locus,
                    group = 1
                ),
                locus = str_extract(.data$a,
                    sample_locus,
                    group = 2
                )
            ) |>
            select(-.data$a) |>
            pivot_wider(
                id_cols = "sample",
                names_from = "locus",
                values_from = "records",
                values_fill = 0
            )
    } else if (inherits(path, "data.frame")) {
        if (!all.variants) {
            message("'var_id' is ignored when 'all.variants' is FALSE.")
            reads <-
                pivot_wider(path,
                    id_cols = "sample",
                    names_from = "locus",
                    values_from = "reads",
                    values_fn = sum,
                    values_fill = 0
                )
        } else if (all.variants) {
            reads <-
                pivot_wider(
                    path,
                    id_cols = "sample",
                    names_from = all_of(c("locus", var_id)),
                    values_from = "reads",
                    values_fill = 0
                )
        }
    }
    return(reads)
}
