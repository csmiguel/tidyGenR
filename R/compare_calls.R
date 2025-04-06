#' Comparision between multiple tidy variants or genotypes.
#'
#' @param td List of named tidy dataframes with at least 'sample',
#' 'locus', 'sequence'.
#' @param dest_file Path to EXCEL file to write results. Default (FALSE),
#' no file is written.
#' @param creads Account for reads. dataframes in 'td' must have 'reads'
#'  column with the
#' number of reads supporting each variant. If FALSE (default)
#' distances between calls
#' and plots are produced on the basis of presence (1)/absense (0) of
#'  the variant.
#' If TRUE, read count are accounted for estimating distances and producing
#'  the plots. Default (FALSE).
#' @returns A list of 4 elements:
#' - *tidy_dat*: combined tidy number of reads for all elements in 'td'. If
#' 'creads = FALSE', cells are coded as '1' (presence) and '0' (absence).
#' - *variants*: list with dataframes in wide format with all variants as
#'  columns.
#'  Variants are named 'locus_\{first 6 characters of md5 of DNA seq\}'.
#'  Zeroes '0'
#'  are included where no variant was found. If
#' 'creads = FALSE', cells are coded as '1' (presence) and '0' (absence)
#' - *comp*: T/F matrix. samples x variants. TRUE, the cell*ij* across all
#' matrices has the same value; FALSE, any of them has a different value.
#' - *dist*: pairwise distances between matrices (analogue to Manhattan
#'  distance).
#' - *plot1*: heatmap with frequency of calls (>1 == 1) across elements
#'  compared.
#' - *plot2*: heatmap with calls for all elements compared.
#'    Black tiles indicate calls absent in the element but present in another.
#' element.
#' - *plot3*: boxplot of reads ordered by frecuency.
#' (output only if 'creads' == TRUE).
#' - *plot4*: MDS plot (output only when length(td) > 2).
#' - (it writes an EXCEL file with the results to 'dest_file'.)
#' @rdname compare_calls
#' @export
#' @examples
#' data("variant_calls")
#' compare_calls(variant_calls, creads = TRUE)
compare_calls <- function(td, dest_file = FALSE, creads = FALSE) {
    # validate input
    cc_validate_td(td)
    # element 1: combined tidy data.
    elem1 <- cc_elem1_combine(td, creads)
    # element 2: list of dataframes: each dataframe for a method where rows
    # are samples and cols locus_variant. Cells are counts.
    elem2 <- cc_elem2_wide_variants(td, elem1)
    # element 3: T/F matrix
    elem3 <- cc_elem3_tf(elem2, creads)
    # element 4: distance between matrices
    elem4 <- cc_elem4_dist(elem2, creads)
    # tidy elem1
    elem1t <-
        elem1 |>
        pivot_longer(
            cols = -c("sample", "locus", "md5"),
            names_to = "method",
            values_to = "value"
        ) |>
        mutate(loc_md5 = glue("{locus}_{sprintf('%.6s', md5)}")) |>
        arrange(.data$sample, .data$loc_md5)
    # plots
    p1 <- cc_p1(elem1)
    p2 <- cc_p2(elem1t, creads)
    p3 <- cc_p3(elem1t, creads)
    p4 <- cc_p4(elem4)
    # write to excel
    cc_write_excel(elem2, elem3, elem4, dest_file)
    # return list
    to_return <-
        list(
            tidy_dat = elem1,
            variants = elem2,
            comp = elem3,
            dist = elem4,
            plot1 = p1,
            plot2 = p2,
            plot3 = p3,
            plot4 = p4
        )
    return(to_return)
}

# helper functions for compare_calls:
#' Validate input for compare calls.
#' @rdname compare_calls
#' @return STOPs if 'td' is not valid for 'compare_calls()'.
#' @noRd
cc_validate_td <- function(td) {
    mand_cols <- c("sample", "locus", "sequence")
    # checks
    stopifnot(
        is.list(td),
        !is.null(names(td)),
        length(td) > 1,
        all(vapply(td, function(x) inherits(x, "data.frame"), logical(1))),
        all(vapply(
            td,
            function(x) {
                all(mand_cols %in% names(x))
            },
            logical(1)
        ))
    )
    message(
        "'td' is a named list with 2 or more dataframes that contain ",
        paste(mand_cols, collapse = " "), "."
    )
}

#' Combine data from multiple tidy variants to compare into one.
#' @rdname compare_calls
#' @return 'dataframe' with 'sample' 'locus' 'md5' and one column per
#' comparison. Each cell under comparisons has 'reads' values. '0' are
#' added to cells in cases with no reads.
#' @noRd
cc_elem1_combine <- function(td, creads) {
    # names list to object
    n1 <- names(td)
    # element 1: combined tidy data.
    # add md5 and select only sample, locus, md5
    vdigest <- Vectorize(digest)
    if (creads) {
        z <-
            n1 |>
            lapply(function(x) {
                td[[x]] |>
                    rowwise() |>
                    mutate(md5 = vdigest(.data$sequence, "md5")) |>
                    group_by(.data$locus, .data$sample, .data$md5) |>
                    summarise(
                        dset = x,
                        reads = .data$reads,
                        .groups = "drop"
                    ) |>
                    ungroup()
            }) |>
            do.call(what = "rbind") |>
            pivot_wider(
                id_cols = c("sample", "locus", "md5"),
                names_from = "dset",
                values_from = "reads",
                values_fn = sum,
                values_fill = 0
            )
    } else if (!creads) {
     z <-
        n1 |>
        lapply(function(x) {
            td[[x]] |>
                rowwise() |>
                mutate(md5 = vdigest(.data$sequence, "md5")) |>
                group_by(.data$locus, .data$sample, .data$md5) |>
                summarise(dset = x, .groups = "drop") |>
                ungroup()
        }) |>
        do.call(what = "rbind") |>
        pivot_wider(
            id_cols = c("sample", "locus", "md5"),
            names_from = "dset",
            values_from = "dset",
            values_fn = length,
            values_fill = 0
        )
    }
    return(z)
}

#' Combine data from multiple tidy variants to compare into one.
#'
#' Used inside 'compare_calls()'.
#' @rdname compare_calls
#' @param elem1 Output from 'cc_elem1_combine()'.
#' @return List of 'dataframes' with 'sample' and one
#' column per 'locus_first_6_char_MD5'.
#' Cells are read count.
#' @noRd
cc_elem2_wide_variants <- function(td, elem1) {
    # names from list
    n1 <- names(td)
    n1 |>
        lapply(function(x) {
            elem1 |>
                select(.data$sample, .data$locus, .data$md5, x) |>
                pivot_wider(
                    id_cols = "sample",
                    names_from = c("locus", "md5"),
                    names_glue = "{locus}_{sprintf('%.6s', md5)}",
                    values_from = x,
                    names_sort = TRUE,
                    values_fill = 0
                )
        }) |>
        setNames(n1)
}

#' Create T/F tablmae
#'
#' Used inside 'compare_calls()'.
#' @rdname compare_calls
#' @param elem2 Output from 'cc_elem2_wide_variants()'.
#' @return 'dataframe' with 'sample' one column per 'locus_first_6_char_MD5'.
#' Cells are TRUE if same call across all methods, or FALSE if any of the calls
#'  differs.
#' @noRd
cc_elem3_tf <- function(elem2, creads) {
    # convert all dataframes in list to matrices
    w <-
        lapply(elem2, function(x) as.matrix(x[, -1]))
    # convert to an array
    y <-
        array(unlist(w), dim = c(
            nrow(w[[1]]),
            ncol(w[[1]]),
            length(w)
        ))
    # true false if cell_j_i equal across array
    z <-
        apply(y, c(1, 2), function(x) all(x == x[1]))
    if (creads) {
        z <-
            apply(y, c(1, 2), function(x) all(x > 0) | all(x == 0))
    }
    colnames(z) <- names(elem2[[1]][-1])
    rownames(z) <- elem2[[1]]$sample
    return(z)
}

#' Create table with pairwise distances
#'
#' Used inside 'compare_calls()'.
#' @rdname compare_calls
#' @return 'dataframe' with 'sample' and one column per
#' 'locus_first_6_char_MD5'.
#' Cells are TRUE if same call across all methods, or FALSE if any of the calls
#'  differs.
#' @noRd
cc_elem4_dist <- function(elem2, creads) {
    z <- dist_m_all(elem2)
    if (creads) {
        w <-
            lapply(elem2, function(x) as.matrix(x[, -1]))
        z <-
            dist_m_all(lapply(w, function(x) {
                apply(x, 2, function(y) y / (sum(y) + 1))
            }))
    }
    return(z)
}

#' Write results from compare_calls() to EXCEL file
#'
#' Used inside 'compare_calls()'.
#' @rdname compare_calls
#' @param elem3 Output from cc_elem3_tf().
#' @param elem4 Output from cc_elem4_dist().
#' @return EXCEL file with sheets having elem2, elem3, elem4.
#' @noRd
cc_write_excel <- function(elem2, elem3, elem4, dest_file) {
    if (isFALSE(dest_file)) {
        message(
            "No output EXCEL file has been written. ",
            "Change this behaviour by setting a path to 'dest_file'."
        )
    } else if (is.character(dest_file)) {
        if (file.exists(dest_file)) {
            file.remove(dest_file)
            }
        zz <-
            ldply(elem2, .id = "method") |>
            as_tibble() |>
            arrange(.data$sample, .data$method)
        l <-
            c(
                list("all_methods" = zz),
                elem2,
                list(
                    "comparisonTF" = as.data.frame(elem3),
                    "distance" = elem4
                )
            )
        write_xlsx(l, dest_file)
    }
}

#' Heatmap with frecuency of each variant call across methods compared
#'
#' Used inside 'compare_calls()'.
#' @rdname compare_calls
#' @return Heatmap.
#' @noRd
cc_p1 <- function(elem1) {
    elem1_temp <-
        elem1 |>
        mutate_if(is.numeric, function(x) as.numeric(x > 0))
    elem1_01 <-
        elem1_temp |>
        mutate(sum01 = rowSums(select(
            elem1_temp,
            where(is.numeric)
        ))) |>
        mutate(loc_md5 = glue("{locus}_{sprintf('%.6s', md5)}")) |>
        arrange(.data$sample, .data$loc_md5)
    z <-
        ggplot(
            elem1_01,
            aes(
                x = .data$sample,
                y = .data$loc_md5,
                fill = as.factor(.data$sum01)
            )
        ) +
        geom_tile(color = "white") +
        scale_fill_brewer(name = "frequency", palette = "Set2") +
        theme_classic() +
        theme(
            axis.text.x = element_text(
                angle = 70,
                vjust = 1,
                hjust = 1,
                size = 6
            ),
            axis.text.y = element_text(size = 6),
            legend.position = "right",
            legend.text = element_text(
                colour = "black",
                size = 8
            )
        )
    return(z)
}

#' Heatmaps comparing methods
#'
#' Used inside 'compare_calls()'. Red tiles are absent calls, but
#' which are present in other methods. Colors according to
#' presence/absence ('creads = FALSE') or read count ('creads = TRUE').
#' @rdname compare_calls
#' @param elem1t Tranformed elem1 (see 'compare_calls()').
#' @return Facetted heatmaps.
#' @noRd
cc_p2 <- function(elem1t, creads) {
    if (creads) {
        z <-
            ggplot(
                elem1t,
                aes(
                    x = .data$sample,
                    y = .data$loc_md5,
                    fill = .data$value
                )
            ) +
            scale_fill_gradientn(
                trans = "log",
                colors = c("red", "orange", "yellow", "darkgreen"),
                name = "Reads (log scale)",
                breaks = c(0, 1, 10, 100, 1000, 10000),
                na.value = "black"
            )
    } else if (!creads) {
        z <-
            ggplot(
                elem1t,
                aes(
                    x = .data$sample,
                    y = .data$loc_md5,
                    fill = as.factor(.data$value)
                )
            ) +
            scale_fill_manual(
                name = "presence",
                values = c(
                    "1" = "green",
                    "0" = "black"
                )
            )
    }
    z <-
        z +
        geom_tile(color = "white") +
        facet_wrap(~ .data$method) +
        theme_classic() +
        theme(
            axis.text.x = element_text(
                angle = 70,
                vjust = 1,
                hjust = 1,
                size = 6
            ),
            axis.text.y = element_text(size = 6),
            legend.position = "right",
            legend.text = element_text(
                colour = "black",
                size = 8
            )
        )
    return(z)
}

#' Boxplots comparing number of reads between methods
#'
#' Used inside 'compare_calls()'.
#' @rdname compare_calls
#' @return Boxplots with loci ordered by read number.
#' @noRd
cc_p3 <- function(elem1t, creads) {
    if (creads) {
        z <-
            elem1t |>
            mutate(locus = reorder(.data$locus, .data$value)) |>
            ggplot(aes(
                x = .data$locus,
                y = .data$value,
                fill = .data$method
            )) +
            geom_boxplot(
                linewidth = 0.3,
                outlier.size = .3,
                color = "grey20"
            ) +
            scale_y_continuous(trans = "log10") +
            ylab("Number of reads") +
            xlab("Locus") +
            theme_classic() +
            theme(
                axis.text.x = element_text(
                    angle = 70,
                    vjust = 1,
                    hjust = 1,
                    size = 6
                ),
                axis.text.y = element_text(size = 6)
            )
    } else if (!creads) {
        z <- NULL
    }
    return(z)
}

#' MDS distance between methods
#'
#' Used inside 'compare_calls()'. Computes pairwise distance
#' between methods and plots MDS. (see 'dist_m' internally).
#' @rdname compare_calls
#' @return MDS plot.
#' @noRd
cc_p4 <- function(elem4) {
    if (nrow(elem4) > 1) {
        z <- mds_comp(elem4)
    } else if (nrow(elem4) == 1) {
        z <- NULL
        message("MDS with 2 dimensions requires at least 3 elements to compare.")
    }
    return(z)
}
