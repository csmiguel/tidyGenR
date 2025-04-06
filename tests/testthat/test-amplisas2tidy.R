test_that("amplisas2tidy works", {
    data("variants")
    fp <-
        list.files(system.file("extdata", "amplisas", package = "tidyGenR"),
            full.names = TRUE
        )
    tidy_vars <- amplisas2tidy(fp[1:3])
    expect_equal(names(variants), names(tidy_vars))
})
