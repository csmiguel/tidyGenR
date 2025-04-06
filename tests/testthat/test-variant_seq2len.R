test_that("var_seq2len returns a data.frame with given columns and variants coded as integers.", {
    data(variants)
    v <- var_seq2len(variants)
    vvariant <- v$variant
    expect_equal(class(vvariant), "integer")
    expect_equal(colnames(v), c("locus", "sample", "variant", "reads"))
    #
})

test_that("var_seq2len fails when malformed dataframe.", {
    data(variants)
    variants <- dplyr::rename(variants, malformed = locus)
    expect_error(var_seq2len(variants))
})
