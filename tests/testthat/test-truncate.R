# Download per-locus demultiplexed fastq
dem <-
  system.file("extdata", "demultiplexed", package = "tidyGenR")

# default parameters work
test_that("Default parameters work", {
    loci <- c("chrna9", "nfkbia")
    tr <-
        trunc_amp(
            loci = loci,
            mode_trun = "pe",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            rv_pattern = "2.fastq.gz",
            trunc_fr = 200,
            max_ee = c(3, 3),
            outdir = file.path(tempdir(), "truncated")
        )
    expect_equal(names(tr), loci)
    expect_equal(class(tr), "list")
    expect_true("matrix" %in% class(tr[[1]]))
    expect_equal(ncol(tr[[1]]), 2)
})

test_that("paired-end accepts arguments correctly. paired-end, all loci:", {
    tr_pe <-
        trunc_amp(
            mode_trun = "pe",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            rv_pattern = "2.fastq.gz",
            trunc_fr = 200,
            max_ee = c(3, 3),
            outdir = file.path(tempdir(), "truncated")
        )
    # pe retturns error when rv_pattern is not stated
    expect_error(
        trunc_amp(
            loci = "chrna9",
            mode_trun = "pe",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            # rv_pattern = "2.fastq.gz",
            trunc_fr = 200,
            max_ee = c(3, 3),
            outdir = file.path(tempdir(), "truncated")
        )
    )
    # max_ee of numeric(1) is equivalent to numeric(2)
    tr_pe_maxee1 <-
        trunc_amp(
            loci = "chrna9",
            mode_trun = "pe",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            rv_pattern = "2.fastq.gz",
            trunc_fr = 200,
            max_ee = 3,
            outdir = file.path(tempdir(), "truncated")
        )
    tr_pe_maxee2 <-
        trunc_amp(
            loci = "chrna9",
            mode_trun = "pe",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            rv_pattern = "2.fastq.gz",
            trunc_fr = 200,
            max_ee = c(3, 3),
            outdir = file.path(tempdir(), "truncated")
        )
    expect_identical(tr_pe_maxee1, tr_pe_maxee2)

    # trunc_fr of numeric(1) is equivalent to numeric(2)
    tr_pe1 <-
        trunc_amp(
            loci = "chrna9",
            mode_trun = "pe",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            rv_pattern = "2.fastq.gz",
            trunc_fr = 200,
            max_ee = 3,
            outdir = file.path(tempdir(), "truncated")
        )
    tr_pe2 <-
        trunc_amp(
            loci = "chrna9",
            mode_trun = "pe",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            rv_pattern = "2.fastq.gz",
            trunc_fr = c(200, 200),
            max_ee = 3,
            outdir = file.path(tempdir(), "truncated")
        )
    expect_identical(tr_pe2, tr_pe1)
})

# out folder is created if it does not exist
test_that("folder in 'outdir' is created if it does not exist", {
    unlink(file.path(tempdir(), "truncated"), recursive = T)
    expect_message(
        trunc_amp(
            loci = "chrna9",
            mode_trun = "se",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            trunc_fr = 200,
            outdir = file.path(tempdir(), "truncated")
        ),
        regexp = ".*Creating output directory.*"
    )
    expect_true(dir.exists(file.path(tempdir(), "truncated")))
})


# only first numeric is taken when fed with numeric(2) trunc_fr
test_that("only first numeric is taken when fed with numeric(2) trunc_fr in single end", {
  expect_no_error(
        trunc_amp(
            loci = "chrna9",
            mode_trun = "se",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            trunc_fr = c(100, 200),
            outdir = file.path(tempdir(), "truncated")
        ))
  }
  )

# returns error when se is fed with numeric(2) in max_ee
test_that("not correct length of max_ee", {
    expect_error(
        trunc_amp(
            loci = "chrna9",
            mode_trun = "se",
            in_dir = dem,
            fw_pattern = "1.fastq.gz",
            trunc_fr = 200,
            max_ee = c(3, 3),
            outdir = file.path(tempdir(), "truncated")
        )
    )
})

# trunc_fr dataframe feeds with locus-specific truncation lengths
test_that("trunc_fr dataframe feeds with locus-specific truncation lengths. trunc_amp() returns a list of length 2", {
  data("trunc_fr")
  h <- trunc_amp(
    loci = c("chrna9", "nfkbia"),
    mode_trun = "pe",
    in_dir = dem,
    fw_pattern = "1.fastq.gz",
    rv_pattern = "2.fastq.gz",
    trunc_fr = trunc_fr,
    outdir = file.path(tempdir(), "truncated")
  )

  expect_length(h, 2)  # Check if the list has exactly 2 elements
})

