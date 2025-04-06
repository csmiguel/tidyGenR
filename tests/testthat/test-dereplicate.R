truncated <- system.file("extdata", "truncated", package = "tidyGenR")
fs <- list.files(truncated,
    pattern = "fastq.gz", full.names = T
)

# default works
test_that("expected list is produced with default parameters", {
    ln <- 5
    dr <-
        dereplicate(fs[1:ln])
    # result is a list
    expect_equal(class(dr), "list")
    # length is as expected
    expect_length(dr, ln)
    # all elements in list are data.frames
    expect_true(all(vapply(dr, function(x) "data.frame" %in% class(x), logical(1))))
    # list names == extracted groups, are OK
    expect_equal(
        names(dr),
        c("chrna9_F", "chrna9_R", "nfkbia_F", "nfkbia_R", "rogdi_F")
    )
})

# min_sam_fr filter (filter cells below a number) work
test_that("'min_sam_fr' filter works:", {
    ln <- 5
    dr_0 <- dereplicate(fs[1:ln], min_sam_fr = 0, min_loc_fr = 0)
    min_sam_fr_100 <- 10
    # the minimum value in a cell is > 10
    dr_100 <-
        dereplicate(fs[1:ln], min_sam_fr = min_sam_fr_100, min_loc_fr = 0)
    expect_true(all(unlist(plyr::llply(
        dr_100,
        function(x) {
            all(x[, 3, drop = TRUE] > min_sam_fr_100)
        }
    ))))
    # any of the values from default is below min_sam_fr_100
    expect_true(all(unlist(plyr::llply(
        dr_0,
        function(x) {
            any(x[, 3, drop = TRUE] < min_sam_fr_100)
        }
    ))))
})

# min_loc_fr filter (filter cells below a number) work
test_that("'min_sam_fr' filter works:", {
    lc <- c("chrna9", "nfkbia")
    ln <- fs[grep(paste(lc, collapse = "|"), fs)]

    dr_0 <- dereplicate(ln,
        min_sam_fr = 0,
        min_loc_fr = 0
    )
    min_loc_fr10 <- 10
    # the minimum rowSums (variant count across all samples) is above
    dr_10 <-
        dereplicate(ln,
            min_sam_fr = 0,
            min_loc_fr = min_loc_fr10
        )
    rw_sums <-
        lapply(dr_10, function(x) {
            z <-
                select(x, -c(1, 2)) |>
                rowSums()
            all(z > min_loc_fr10)
        })
    expect_true(all(unlist(rw_sums)))
})

# 'by' argument selects field correctly
test_that("'by' for selecting another field works", {
    lc <- c("chrna9", "nfkbia")
    ln <- fs[grep(paste(lc, collapse = "|"), fs)]
    # by sample
    dr_sample <- dereplicate(ln,
        by = "([a-zA-Z0-9]*)_[a-zA-Z0-9]*_[F|R]"
    )
    expect_equal(names(dr_sample), c("BOR1061", "BOR1063", "BOR1069"))
    # by locus
    dr_locus <- dereplicate(ln,
        by = "[a-zA-Z0-9]*_([a-zA-Z0-9]*)_[F|R]"
    )
    expect_equal(names(dr_locus), lc)
})

# out_xlsx writes excel
test_that("'by' for selecting another field works", {
    # by sample
    suppressWarnings(
        file.remove(file.path(tempdir(), "test.xlsx"))
    )
    dr_excel <- dereplicate(fs[1:5],
        out_xlsx = file.path(tempdir(), "test.xlsx")
    )
    expect_true(file.exists(file.path(tempdir(), "test.xlsx")))
})
