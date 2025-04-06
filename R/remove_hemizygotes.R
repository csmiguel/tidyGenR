#' Remove hemizygote calls
#'
#' Removes all rows with hemizygous alleles from tidy genotypes.
#' @details If "A/B" -> "A/B"; if "A" -> NULL.
#' It only is applicable to diploid data.
#' @param gen Tidy diploid genotypes.
#' @return Tidy genotypes without hemizygous alleles.
#' @examples
#' data("genotypes")
#' remove_hemizygotes(genotypes)
#' @export
remove_hemizygotes <- function(gen) {
    if (guess_ploidy(gen) != 2) {
        stop("Genotypes have to be diploid.")
    }
    z <-
        ddply(gen, ~ sample + locus, function(x) {
            if (is.na(nrow(x))) {
                NULL
            } else if (nrow(x) == 1) {
                NULL
            } else if (nrow(x) == 2) {
                x
            }
        })
    return(as_tibble(z))
}
