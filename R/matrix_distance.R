#'  Distance between two matrices
#'
#' Computed as the sum of absolute difference across cells.
#' Analogue to Manhattan distance.
#'
#' @param m1 Matrix 1.
#' @param m2 Matrix 2.
#' @param prop If TRUE, it divides number of differences by length(matrix).
#' @returns Numeric distance, analogue to Manhattan distance.
dist_m <- function(m1, m2, prop = FALSE) {
    if (!identical(dim(m1), dim(m2))) {
        stop("Dimensions of matrices m1 and m2 must be identical.")
    }
    distm <- sum(abs(m1 - m2))
    if (prop) {
        distm <- distm / length(m1)
    }
    return(distm)
}

#' Applies distance between all tables in a list.
#'
#' Requires a list of matrices in which the first column are samples,
#' and col2...coln are numeric. Ex. number of observed variants in each cell.
#' @param ld List of dataframes.
#' @returns Dataframe with pairwise distances between all combinations.
#' Absolute distance 'dist_euc' and relative 'dist_eucp'.
dist_m_all <- function(ld) {
    # reference df with pairwise comparisons.
    e4 <-
        as.data.frame(t(combn(names(ld), 2)))
    names(e4) <- c("method1", "method2")
    # for each pairwise comparison compute distance.
    z <-
        apply(e4, 1, function(x) {
            m1 <- as.matrix(ld[[as.character(x["method1"])]][, -1])
            m2 <- as.matrix(ld[[as.character(x["method2"])]][, -1])
            de <- dist_m(m1, m2)
            dep <- dist_m(m1, m2, prop = TRUE)
            return(c(de, dep))
        })
    z <- data.frame(t(z))
    names(z) <- c("dist_euc", "dist_eucp")
    return(cbind(e4, z))
}
