truncated <-
  system.file("extdata", "truncated", package = "tidyGenR")

test_that("reads are tracked", {
    # data frame with fake variants
    data("variants")
    # named list with elements to track
    l <-
        list(
            truncated = c(truncated, "F_filt.fastq.gz"),
            variants = variants
            )
    rt <- reads_track(l)
    # names and class are ok
    expect_equal(names(rt), c("sample", names(l)))
    expect_true(any(class(rt) %in% "data.frame"))
    # list of size 1 works
    l1 <-
        list(
            truncated = c(truncated, "F_filt.fastq.gz")
        )
    expect_no_error(reads_track(l))
})
