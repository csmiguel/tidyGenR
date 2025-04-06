test_that("returns a list of equal length as input", {
    data(qlist)
    l1 <- length(qlist)
    expect_equal(l1, length(tidyGenR:::align_matrices(qlist)))
    expect_equal(class(tidyGenR:::align_matrices(qlist)), "list")
})

test_that("returns error when leap between consecutive ordered Ks > 1", {
    data(qlist)
    qlist_cut <- qlist[vapply(qlist, ncol, numeric(1)) != 2]
    expect_error(align_matrices(qlist_cut),
        "The maximum step increment in no cols in 1.",
        fixed = TRUE
    )
})
