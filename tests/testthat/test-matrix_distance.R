test_that("already tested in compare_calls()", {
    expect_true(
        exists("compare_calls", mode = "function", where = "package:tidyGenR")
    )
})
