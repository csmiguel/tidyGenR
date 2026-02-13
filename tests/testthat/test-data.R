test_that("genotypes is formatted correctly", {
    data("genotypes")
    expect_true(
        all(
            c("sample", "locus", "allele", "allele_no", "reads", "sequence") %in%
                names(genotypes)
        )
    )
})

test_that("variants is formatted correctly", {
    data("variants")
    expect_true(
        all(
            c("sample", "locus", "variant", "reads", "sequence") %in%
                names(variants)
        )
    )
})

test_that("variant_calls is a list of dataframes formatted correctly", {
    data("variant_calls")
    expect_equal(class(variant_calls), "list")
    expect_true(
        all(vapply(variant_calls, function(x) "data.frame" %in% class(x), logical(1)))
    )
})

test_that("trunc_fr is a dataframe with truncation lengths", {
    data("trunc_fr")
    expect_equal(
        names(trunc_fr),
        c("locus", "trunc_f", "trunc_r")
    )
})


test_that("primers are locus specific primerss", {
    data("primers")
    expect_equal(
        names(primers),
        c("locus", "fw", "rv")
    )
})
