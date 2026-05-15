test_that("default works.", {
    data("variant_calls")
    cc <-
        compare_calls(variant_calls, creads = TRUE)
    # names of ouput list are ok
    expect_equal(names(cc), c(
        "tidy_dat", "variants", "comp", "dist", "plot1",
        "plot2", "plot3", "plot4"
    ))
    # classes of elements in list are ok
    expect_true("data.frame" %in% class(cc[[1]]))
    expect_true("list" %in% class(cc[[2]]))
    expect_true("matrix" %in% class(cc[[3]]))
    expect_true("data.frame" %in% class(cc[[4]]))
    expect_equal(
        unique(vapply(cc[5:length(cc)], function(x) class(x)[2], character(1))),
        "ggplot"
    )
})
# deeper into each object
test_that("each element returns as expected", {
    data("variant_calls")
    cc <-
        compare_calls(variant_calls, creads = TRUE)
    # element 1 : classes of elements in list are ok
    expect_equal(names(cc[[1]]), c(
        "sample", "locus", "md5",
        names(variant_calls)
    ))
    # element 2:
    expect_equal(names(cc[[2]]), names(variant_calls))
    expect_true(all(vapply(
        cc[[2]], function(x) "data.frame" %in% class(x),
        logical(1)
    )))
    # element 3:
    expect_true(is.logical(as.vector(cc[[3]])))
    # element 4:
    expect_equal(names(cc[[4]]), c("method1", "method2", "dist_Manh", "dist_Manhp"))
})

# what about when there is no 'reads'
test_that("when no 'reads' field, calculations are done as expected:", {
    # distance matrix is affected
    data("variant_calls")
    cc <-
        compare_calls(variant_calls, creads = FALSE)
    # only check elements which change respect to TRUE.
    # element1: there are only 0 and 1 in element 1. Indicating 1 when present
    # and 0 absence.
    only_01 <-
        cc[[1]][, 4:(length(variant_calls) + 3)] |>
        unlist() |>
        as.numeric() |>
        unique()
    expect_equal(c(1, 0), only_01)
    # element 2. there are only 0 and 1 in specific dataframes.
    elem2 <-
        vapply(cc[[2]], function(x) {
            z <-
                unlist(x[, -1]) |>
                as.numeric() |>
                unique() |>
                sum()
            z == 1
        }, logical(1))
    expect_true(all(elem2))
})

# excel file is written

test_that("excel file is written", {
    ex_path <- file.path(tempdir(), "cc.xlsx")
    cc <-
        compare_calls(variant_calls, dest_file = ex_path)
    expect_true(file.exists(ex_path))
})

# trigger errors
test_that("excel file is written", {
    # requires more than 2 df
    vars2 <- variant_calls[1]
    expect_error(compare_calls(vars2))
    # requires reads when 'creads = T'
    vars_no_reads <-
        variant_calls |>
        lapply(function(x) select(x, -reads))
    expect_error(
        compare_calls(vars_no_reads[c(1:3)], creads = TRUE)
    )
})
