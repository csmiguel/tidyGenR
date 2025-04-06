#' Remove monorphic loci
#'
#' Removes monomorphic (or polymorphic) loci from tidy genotypes.
#'
#' @param gen Tidy genotypes (output from 'genotype').
#' @param invert T/F, if T, polymorphic loci are removed.
#' @return Tidy genotypes with mono/poly-morphic loci.
#' @examples
#' data(genotypes)
#' remove_monomorphic(genotypes)
#' @export
remove_monomorphic <- function(gen, invert = FALSE) {
    mono <-
        dlply(
            gen, ~locus,
            function(x) length(unique(x$allele))
        ) == 1
    loci_string <- paste(names(mono)[mono], collapse = " ")
    message("Loci ", loci_string, " are monomorphic.")
    if (invert) {
        z <- gen[gen$locus %in% names(mono)[mono], ]
        n_mono <- length(unique(z$locus))
        message("A total of ", n_mono, " monomorphic loci have been returned.")
    } else if (!invert) {
        z <- gen[!gen$locus %in% names(mono)[mono], ]
        n_poly <- length(unique(z$locus))
        message("A total of ", n_poly, " polymorphic loci have been returned.")
    }
    return(z)
}
