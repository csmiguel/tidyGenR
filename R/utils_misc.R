#' Align and view DNA sequences
#' Align and view DNA sequences based on *DECIPHER* functions.
#' @param seqs Character vector with DNA sequences.
#' @param browse Display msa in web browser.
#' @param out_bs Return biostring MSA.
#' @param ... other parameters passed to 'BrowseSeqs()'.
#' @return  Biostring object with aligned sequences and display msa in web
#'  browser.
#' @examples
#' tidyGenR:::view_msa(c("CAAAA", "CAAAT"))
view_msa <- function(seqs, browse = TRUE, out_bs = TRUE, ...) {
    s <- DNAStringSet(seqs)
    sal <- AlignSeqs(s)
    if (browse) {
        BrowseSeqs(sal, ...)
    }
    if (out_bs) {
        return(sal)
    }
}

#' Align Q matrices
#'
#' Align list of q-matrices (matrices with values from 0 to 1) with consistent
#' row numbers. It is meant for the alignment q-matrices derived from
#' STRUCTURE (Pritchard et al. 2000) or similar software.
#' @details
#' It can handle matrices with similar and different number of columns.
#' Matrices in the list are sorted by increasing number of columns and
#' the clusters from the first matrix are used as a reference for
#' ordering clusters from successive q-matrices
#' in the list.
#' @references
#'  Pritchard, J. K., Stephens, M., & Donnelly, P. (2000).
#'  _Inference of population structure using multilocus genotype data_
#'  Genetics, 155(2), 945-959.
#' @param mat List of dataframes or matrices with q-values.
#' @return List of dataframes with aligned clusters.
#' @examples
#' data("qlist")
#' tidyGenR:::align_matrices(qlist)
align_matrices <- function(mat) {
    # mat, list of matrices
    mat <- lapply(mat, as.matrix)
    # order by K
    qk <- mat[order(vapply(mat, ncol, integer(1)))]
    # confirm max step increment is 1
    kn <- vapply(qk, ncol, integer(1)) # no of cols for each matrix
    knc <- as.numeric(rep(NA, length(kn) - 1))
    # increments in ncol from consecutive matrices in the list
    for (i in 2:length(kn)) {
        knc[i - 1] <- kn[i] - kn[i - 1]
    }
    if (max(knc) > 1) {
        stop("The maximum step increment in no cols in 1.")
    }

    ord_matrices <- qk[1]

    for (j in seq_along(qk)[-1]) {
        targetk <- qk[[j]]
        refk <- ord_matrices[[j - 1]]
        no_cols_ref <- ncol(refk)
        no_cols_target <- ncol(targetk)
        col_ord <- rep(NA, no_cols_target)
        for (i in seq_len(no_cols_ref)) {
            vec <- refk[, i]
            col_dists <- apply(targetk, 2, function(x) sum((x - vec)^2))
            best_dist <- min(col_dists[is.na(col_ord)])
            best_col <- match(best_dist, col_dists)[1]
            col_ord[best_col] <- colnames(refk)[i]
        }
        col_ord[is.na(col_ord)] <- paste0("Cluster", no_cols_target)
        colnames(targetk) <- col_ord
        # order columns according to clusters
        targetk <- targetk[, order(colnames(targetk)), drop = FALSE]
        ord_matrices[[j]] <- targetk
    }
    names(ord_matrices) <- names(qk)
    return(lapply(ord_matrices, as.data.frame))
}


#' Filter Multiple Sequence Alignment
#'
#' Selects columns (MSA positions) and rows (DNA sequences) in
#' "DNAMultipleAlignment" or "DNAStringSet" objects.
#' Default values will return all positions and columns.
#'
#' @param msa Multiple sequence alignment stored as a 'Biostrings'
#'  "DNAMultipleAlignment" or "DNAStringSet".
#' @param rows Numeric vector for row (sequence) selection.
#' @param cols Numeric vector for cols (columns in MSA) selection.
#' @examples
#' data("msa")
#' tidyGenR:::filt_msa(msa, 1:2, -c(1:100))
filt_msa <- function(msa, rows = seq_along(msa), cols = NULL) {
    allowed_class <- c("DNAMultipleAlignment", "DNAStringSet")
    if (!class(msa) %in% allowed_class) {
        stop("'msa' must be a ", allowed_class)
    }
    if (is.null(cols)) {
        if (class(msa) == allowed_class[1]) {
            ncols <- seq_len(ncols(msa))
        } else if (class(msa) == allowed_class[2]) {
            wmsa <- unique(width(msa))
            if (length(wmsa) > 1) {
                stop("Multiple sequence lengths in MSA.")
            }
            cols <- seq_len(wmsa)
        }
    }
    alm <- as.matrix(msa)
    almfilt <- alm[rows, cols, drop = FALSE]
    al_filt <-
        DNAStringSet(apply(almfilt, 1, paste0, collapse = ""))
    return(al_filt)
}

#' Remove Ns introduced in 'mergePairs(justConcatenate = T)'
#'
#' Detect N strings and remove them from aligned BStringSet MSA.
#'
#' @details When using 'mergePairs(justConcatenate = T)', a string of 10
#' Ns is introduced between paired-end concatenated reads.
#' This function detects the coordinates of the N string in the MSA and removes
#' it. It uses 'matchPattern()' to detect patterns.
#' WARNING: it has not been tested for multiple hits of the pattern.
#' It could potentially work with other patterns in addition to N strings,
#' but it has not been tested.
#' @returns Filtered BStringSet MSA.
#' @param msa Multiple sequence alignment stored as a "DNAMultipleAlignment"
#'  or "DNAStringSet".
#' @param pattern_n N string to match.
#' @param outfasta Character vector with path to write FASTA file with MSA.
rm_n_msa <- function(msa,
                     pattern_n = paste0(rep("N", 10), collapse = ""),
                     outfasta = FALSE) {
    # detect patterns
    n_ranges <-
        lapply(msa, function(x) {
            matchPattern(pattern_n, x)
        })

    # get dataframe with starting position of N's and their width
    detect_n <-
        seq_along(names(n_ranges)) |>
        ldply(function(x) {
            z <- n_ranges[[x]]@ranges@start
            if (length(z) == 0) {
                z <- "NA"
            }
            w <- n_ranges[[x]]@ranges@width
            if (length(w) == 0) {
                w <- "NA"
            }
            data.frame(
                seqname = names(n_ranges)[x],
                start_n = z,
                width_n = w
            )
        }) |>
        filter(.data$start_n != "NA") |>
        mutate_at(c("start_n", "width_n"), as.numeric, na.rm = TRUE)
    if (nrow(detect_n) == 0) {
        message("\nNo columns have been removed.")
        msa_filt <- msa
    } else if (nrow(detect_n) > 0) {
        # start at
        f <- min(detect_n$start_n)
        # end at
        l <- max(detect_n$start_n + max(detect_n$width_n)) - 1
        # filtered msa
        msa_filt <- filt_msa(msa, cols = -c(f:l))
        message("\nColumns ", f, " to ", l, " have been removed.")
    }
    if (is.character(outfasta)) {
        writeXStringSet(msa_filt, filepath = outfasta)
    }
    return(msa_filt)
}
