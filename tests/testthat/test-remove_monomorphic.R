test_that("there are monomorphic loci and all returned loci are polymorphic", {
    data("genotypes")
    gen_poly <- remove_monomorphic(genotypes)
    nal <-
        select(gen_poly, locus, md5) |>
        distinct() |>
        pull(locus) |>
        table() |>
        as.numeric() |>
        min()
    # min number of alleles in filtered df is > 1
    expect_gt(nal, 1)
    # min number of alleles in non-filtered df is 1
    nal_unfilt <-
        select(genotypes, locus, md5) |>
        distinct() |>
        pull(locus) |>
        table() |>
        as.numeric() |>
        min()
    expect_equal(nal_unfilt, 1)
    # when inverted, all loci are monomorphic
    gen_mono <- remove_monomorphic(genotypes, invert = TRUE)
    nal_mono <-
        select(gen_mono, locus, md5) |>
        distinct() |>
        pull(locus) |>
        table() |>
        as.numeric()
    expect_true(all(nal_mono == 1))
})
