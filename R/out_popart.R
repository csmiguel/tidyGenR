#' Format MSA and traits for PopArt
#'
#' It takes a multiple sequence alignment and data to create a NEXUS file
#'  which can be read by PopART (https://popart.maths.otago.ac.nz).
#' @param msa Multiple sequence alignment, DNAStringSet or DNAMultipleAlignment.
#' @param x Dataframe with trait data.
#' @param blocks NEXUS blocks to add to file. Either "DATA", "TRAITS"
#' or c("DATA", "TRAITS").
#' @param sname Column name in 'x' with sequence names which must be
#' identical to sequence names in the MSA.
#' @param xgroups Column name in 'x' used for creating groups in TRAITS.
#' @param outnex Path to write NEXUS file.
#' @examples
#' data("genotypes")
#' x <-
#'     dplyr::filter(genotypes, locus == "abcg8")
#' nr <- seq_len(nrow(x))
#' y <-
#'     dplyr::mutate(x, sname = paste0("seq_", nr))
#' msa <-
#'     tidy2sequences(y, fasta_header = "{sname}") |>
#'     DECIPHER::AlignSeqs()
#' outnex <- tempfile(fileext = ".nex")
#' out_popart(msa, y,
#'     outnex = outnex,
#'     sname = "sname",
#'     xgroups = "sample",
#'     blocks = c("DATA", "TRAITS")
#' )
#' @export
out_popart <- function(msa,
                       x,
                       blocks = c("DATA", "TRAITS"),
                       sname,
                       xgroups,
                       outnex) {
    # nexus file
    if (file.exists(outnex)) {
        file.remove(outnex)
    }
    #
    if ("DATA" %in% blocks) {
        # class has to be DNAStringSet
        allowedclass <- c("DNAMultipleAlignment", "DNAStringSet")
        if (!inherits(msa, allowedclass)) {
            stop(paste(
                "'msa' must be a multiple sequence alignment in ",
                paste(allowedclass, collapse = " "),
                "format."
            ))
        }
        if (inherits(msa, "DNAMultipleAlignment")) {
            msa <- DNAStringSet(msa)
        }
        # seqs as DNAbin
        seqs_bin <- ape::as.DNAbin(msa)
    }
    if ("TRAITS" %in% blocks) {
        if (!identical(x[[sname]], names(msa))) {
            stop("Names for sequence names and rows in TRAIT  must be identical.")
        }
        # group levels
        tl <-
            unique(x[[xgroups]])
        # trait df
        zz <-
            data.frame(
                pivot_wider(x,
                    id_cols = all_of(sname),
                    names_from = any_of(xgroups),
                    values_fn = length,
                    values_fill = 0,
                    values_from = all_of(xgroups)
                )
            )
        trait_df <-
            unite(zz, "traits",
                all_of(tl),
                sep = ",", remove = TRUE
            )
    }
    # functions to generate blocks
    # data
    fdata <- function() {
        # write DATA block
        ape::write.nexus.data(seqs_bin,
            file = outnex,
            interleaved = FALSE
        )
    }
    # traits
    ftraits <- function() {
        # write TRAITS block
        writeLines(
            c(
                "BEGIN TRAITS;",
                paste0("Dimensions NTRAITS=", length(tl), ";"),
                "Format labels=yes missing=? separator=Comma;",
                paste0("TraitLabels ", paste(tl, collapse = " "), ";"),
                "Matrix"
            ),
            con = CON,
            sep = "\n"
        )
        # append trait
        capture.output(
            write_delim(trait_df,
                file = stdout(),
                delim = " ",
                append = TRUE,
                col_names = FALSE,
                quote = "none",
                progress = FALSE
            ),
            file = CON
        )
        writeLines(c(";", "END;"),
            sep = "\n",
            con = CON
        )
    }
    # cases
    if (all(c("DATA", "TRAITS") %in% blocks)) {
        CON <- file(outnex, "a")
        fdata()
        ftraits()
        close(CON)
    } else if (blocks == "DATA") {
        fdata()
    } else if (blocks == "TRAITS") {
        CON <- file(outnex, "a")
        ftraits()
        close(CON)
    } else {
        stop("No blocks defined.")
    }
    message("NEXUS file written to ", outnex)
}
