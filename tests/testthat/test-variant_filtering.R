test_that("Allelic depth filter and MAF filter work:", {
    data("variants")
    n_unfiltered <- nrow(filter_variants(variants))
    n_filtered <- nrow(filter_variants(variants, ad = 100, maf = 0.1))
    expect_identical(n_unfiltered, nrow(variants))
    expect_identical(n_filtered, 1321L)
})

test_that("output variants have the same columns as input", {
    data("variants")
    fvar <- filter_variants(variants)
    expect_equal(colnames(variants), colnames(fvar))
})

test_that("'invert' returns excluded variants", {
    data("variants")
    n_excluded <- nrow(filter_variants(variants,
        ad = 100,
        maf = 0.1, invert = TRUE
    ))
    n_retained <- nrow(filter_variants(variants,
        ad = 100,
        maf = 0.1, invert = FALSE
    ))
    expect_equal(n_excluded + n_retained, nrow(variants))
})
