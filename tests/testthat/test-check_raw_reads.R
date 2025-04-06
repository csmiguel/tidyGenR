raw <- system.file("extdata", "raw", package = "tidyGenR")
freads <-
    list.files(raw,
        pattern = "1.fastq",
        full.names = TRUE
    )
rreads <-
    list.files(raw,
        pattern = "2.fastq",
        full.names = TRUE
    )
# default single-end and paired-end mode work fine with default parameters.
test_that("se and pe return list with default parameters. All checks passed.", {
    pe <- check_raw_reads(freads = freads, rreads = rreads)
    se <- check_raw_reads(freads = freads)
    expect_equal(vapply(list(pe, se), class, character(1)), c("list", "list"))
    expect_true(all(c(pe$checks, se$checks)))
    expect_equal(pe$samples, c("BOR1061", "BOR1063", "BOR1069"))
})

# returns error and no messages with malformed paths
test_that("returns error when path to fastq are empty or malformed.", {
    expect_error(check_raw_reads(freads = "malformed_path"))
    expect_error(check_raw_reads(
        freads = freads,
        rreads = "malformed_path"
    ))
    expect_error(check_raw_reads())
})

# only check 1 fails
test_that("1 or more files below read threshold", {
    testx <- function() {
        se <-
            check_raw_reads(
                freads = freads,
                low_readcount = 70000
            )
        return(se)
    }

    suppressWarnings(
        expect_false(testx()$checks[1])
    )
    expect_warning(testx())
})

# only check 2 fails.
test_that("error when duplicated sample names", {
    fr <- freads[c(1, 1, 2)]
    expect_error(
        check_raw_reads(freads = fr)
    )
})

# only check 3 fails
test_that("orphan files", {
    fr <- freads[c(1, 2)]
    expect_warning(
        check_raw_reads(
            freads = fr,
            rreads = rreads
        ),
        regexp = "NOT all F files have their corresponding R file."
    )
})

# fails when freads and rreads are both null
test_that("fails when freads and rreads are both NULL", {
    fr <- NULL
    rv <- NULL
    expect_error(
        check_raw_reads(
            freads = fr,
            rreads = rv
        )
    )
})

