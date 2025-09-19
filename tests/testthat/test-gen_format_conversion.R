data("genotypes")
# gen_tidy2compact
test_that("gen_tidy2compact works fine", {
    compact <- gen_tidy2compact(genotypes)
    # col names are as expected
    expect_equal(names(compact), c("sample", "locus", "genotype"))
    # default delim is /.
    expect_true(grepl("/", compact[1, 3]))
    # returns error when tidy has missing one of the needed variables
    gen_wrong_names <- rename(genotypes, anyname = locus)
    expect_error(gen_tidy2compact(gen_wrong_names))
    # any other delim works
    compact2 <- gen_tidy2compact(genotypes, delim = "xxx")
    expect_true(grepl("xxx", compact2[1, 3]))
    # only sample, locus, allele are needed to run
    gen2 <- select(genotypes, sample, locus, allele)
    compact3 <- gen_tidy2compact(gen2)
    expect_equal(compact, compact3) # same as default
})

# gen_compact2wide
test_that("gen_compact2wide works fine", {
    cgen <- gen_tidy2compact(genotypes)
    wgen <- gen_compact2wide(cgen)
    nsamples <- length(unique(genotypes$sample))
    nloci <- length(unique(genotypes$locus))
    # number of rows match no samples
    expect_equal(dim(wgen)[1], nsamples)
    # no of cols match no loci
    expect_equal(dim(wgen)[2] - 1, nloci)
})

# gen_tidy2wide
test_that("gen_tidy2wide works fine", {
    wgen <- gen_tidy2wide(genotypes)
    nsamples <- length(unique(genotypes$sample))
    nloci <- length(unique(genotypes$locus))
    # number of rows match no samples
    expect_equal(dim(wgen)[1], nsamples)
    # no of cols match no loci
    expect_equal(dim(wgen)[2] - 1, nloci)
})

# gen_wide2structure
test_that("gen_wide2structure no error", {
    str1 <- file.path(tempdir(), "test1.str")
    meta <-
        read.csv(system.file("extdata/metadata.csv",
            package = "tidyGenR"
        ))
    # from wide to str
    gen_wide2structure(gen_tidy2wide(genotypes),
        write_out = str1,
        popdata = meta,
        delim = "/"
    )
    expect_snapshot_file(str1, "w2str_snap")
})

# gen_tidy2genalex
test_that("genalex file is formatted correctly:", {
    glx_txt <- file.path(tempdir(), "out.txt")
    glx_xlsx <- file.path(tempdir(), "out.xlsx")
    meta <-
        read.csv(system.file("extdata/metadata.csv",
            package = "tidyGenR"
        ))
    expect_no_error(
        gen_tidy2genalex(genotypes,
            popdata = meta,
            write_out = glx_txt
        )
    )
    expect_no_error(gen_tidy2genalex(genotypes,
        popdata = meta,
        write_out = glx_xlsx
    ))
    expect_snapshot_file(glx_txt, "snap_genalex.txt")
})
# gen_tidy2integers
test_that("alleles are recoded as integers", {
    gen_int <-
        gen_tidy2integers(genotypes)
    expect_true(is.integer(gen_int$allele))
    generr <- rename(genotypes, whatever = allele)
    expect_error(gen_tidy2integers(generr))
})
