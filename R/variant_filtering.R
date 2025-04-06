#' Filter variants based on frequency and depth
#'
#' @param tidy_var Dataframe or tibble with tidy variants.
#' @param maf Miminum variant frequency. For each sample, variants with a
#'  proportion of number of reads across all variants > maf are retained.
#' @param ad  Allele depth. Minimum number of reads supporting a variant.
#' Variants with ad > 3 are retained.
#' @param invert FALSE (default) returns variants passing filters; TRUE inverts
#' the selection.
#' @returns Tidy tibble with filtered variants.
#' @examples
#' data("variants")
#' filter_variants(variants, ad = 100, maf = 0.2)
#' @export
#' @return dataframe with filtered variants/alleles. Variants are named so
#' column names in df are "locus", "sample", "sequence", "reads", "allele".
filter_variants <- function(tidy_var, maf = 0, ad = 0, invert = FALSE) {
    if (!"data.frame" %in% class(tidy_var)) {
        stop("'tidy_var' must be a 'dataframe'.")
    }
    if (!all(c("sample", "locus", "reads") %in% names(tidy_var))) {
        stop("'tidy_var' must have 'sample' 'locus' 'reads'.")
    }
    fvar <-
        ddply(tidy_var, ~ locus + sample,
            .drop = TRUE,
            function(x) {
                z <- x[["reads"]]
                rel <- z / sum(z)
                vmaf <- as.logical(rel > maf)
                vad <- as.logical(z > ad)
                if (!invert) {
                    zz <- x[as.logical(vmaf * vad), ]
                } else if (invert) {
                    zz <- x[!as.logical(vmaf * vad), ]
                }
                return(zz)
            }
        )
    fvar <- as_tibble(drop_na(fvar))
    attr(fvar, "maf") <- maf
    attr(fvar, "ad") <- ad
    return(fvar)
}
