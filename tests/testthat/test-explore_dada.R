truncated <-
  system.file("extdata", "truncated", package = "tidyGenR")
paths <-
    list.files(truncated,
        pattern = "_F_", full.names = TRUE
    )

test_that("default works", {
    set.seed(123)
    fs <- sample(paths, 5)
    ed <-
        explore_dada(fs,
            value_na = 10,
            reduced = TRUE,
            pool = FALSE,
            vline = 2,
            hline_fr = 0.1,
            omega_a = 0.9,
            band_size = 16
        )
    # length is ok
    expect_equal(length(ed), 5L)
    # classes are ok
    expect_true(inherits(ed[[1]], "data.frame"))
    expect_true(all(sapply(2:5, function(i) inherits(ed[[i]], "ggplot"))))
    # tidy dada names are ok
    expect_equal(names(ed$tidy_dada), c(
        "sample", "locus", "sequence",
        "rank", "loglogp", "abundance", "vfreq", "rank2"
    ))
})
