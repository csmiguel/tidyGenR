#' Read AmpliSAT results into tidy variants
#'
#' Reads AmpliSAT results as tidy variants.
#'
#' @details Converts text files in folder with AmpliSAS (Sebastian et al. 2016)
#'  results into a tidy variants dataframe. MD5 is regenerated as the one in
#'  AmpliSAT does not coincide with the MD5 associated with the DNA sequence.
#' @param fp Character vector with paths to 'txt' files from AmpliSAT results.
#' E.x. from 'allseqs', 'filtered', 'clustered'.
#' @references
#' Sebastian et al. (2016). _AMPLISAS: a web server for multilocus genotyping
#'  using next‐generation amplicon sequencing data_.
#' Molecular Ecology Resources, 16(2), 498-510.
#' @return Tidy variants (see 'variant_call()'').
#' @examples
#' fp <- list.files(system.file("extdata", "amplisas", package = "tidyGenR"),
#'     full.names = TRUE
#' )
#' amplisas2tidy(fp)
#' @export
amplisas2tidy <- function(fp) {
    # loci names
    loci <-
        basename(fp) |>
        str_remove(".txt")
    message("\nLoci detected: '", paste(loci, collapse = "', '"), "'.\n")
    message("Files detected: \n", paste(fp, collapse = "\n"))
    vdigest <- Vectorize(digest)
    # apply to all loci: read '.txt' and convert to tidy variants
    z <-
        fp |> # for each locus
        lapply(function(fpx) {
            # read amplisas output
            read.table(fpx,
                sep = "\t",
                header = TRUE
            ) |>
                # reformat
                select(-.data$DEPTH, -.data$SAMPLES, -.data$MD5) |>
                rename(
                    variant = .data$NAME,
                    sequence = .data$SEQ,
                    nt = .data$LENGTH
                ) |>
                rowwise() |>
                mutate(
                    locus = str_split(.data$variant, "-")[[1]][1],
                    variant = str_split(.data$variant, "-")[[1]][2],
                    md5 = vdigest(.data$sequence)
                ) |>
                pivot_longer(
                    cols = -c(
                        "variant", "sequence", "md5",
                        "nt", "locus"
                    ),
                    names_to = "sample",
                    values_to = "reads"
                ) |>
                drop_na() |>
                select(
                    .data$sample, .data$locus, .data$variant, .data$reads,
                    .data$nt, .data$md5, .data$sequence
                )
        }) |>
        do.call(what = "rbind") |>
        arrange(.data$sample, .data$locus, .data$variant)
    return(z)
}
