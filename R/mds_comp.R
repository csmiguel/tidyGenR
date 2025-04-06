#' plot MDS
#'
#' @description
#' Plot MDS from distances between compared calls
#' @param df Dataframe with pairwise distances from 'dist_m()'.
#' One row per unique comparison.
mds_comp <- function(df) {
    if (!is.data.frame(df)) {
        stop("'df' has to bee a dataframe.")
    }
    mand_names <- c("method1", "method2", "dist_euc")
    if (!all(mand_names %in% names(df))) {
        stop("'df' names need to have ", mand_names)
    }
    df <- select(df, all_of(mand_names))
    # Get all unique method names
    all_methods <- unique(c(df$method1, df$method2))
    # Create a complete dataframe with all combinations of method1 and method2
    complete_df <-
        expand.grid(
            method1 = all_methods,
            method2 = all_methods
        ) |>
        left_join(df, by = c("method1", "method2"))
    # Spread to wide format to create the distance matrix
    dist_matrix <-
        complete_df |>
        pivot_wider(
            names_from = "method2",
            values_from = "dist_euc"
        ) |>
        column_to_rownames("method1") |>
        as.matrix()
    # Convert the matrix to a distance object
    dist_object <- as.dist(t(dist_matrix),
        diag = TRUE,
        upper = TRUE
    )
    # Perform Classical MDS
    mds_result <- cmdscale(dist_object, k = 2) # k = 2 for 2D plot
    # Convert result to data frame for plotting
    mds_df <-
        as.data.frame(mds_result) |>
        setNames(c("Dimension1", "Dimension2")) |>
        rownames_to_column("method")
    # Plot the MDS results
    p1 <-
        ggplot(
            mds_df,
            aes(
                x = .data$Dimension1,
                y = .data$Dimension2,
                label = .data$method
            )
        ) +
        geom_point() +
        theme_minimal() +
        labs(
            title = "MDS Plot of Distance Between Matrices",
            x = "Dimension 1",
            y = "Dimension 2"
        )
    return(p1)
}
