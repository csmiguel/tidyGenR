data("variants")
var_red <- variants[c(1:100), ]

# default works correctly
test_that("genotype works with default parameters", {
    gen <-
        genotype(var_red)
    # is a data.frame
    expect_true("data.frame" %in% class(gen))
    # colnames is ok
    expect_true(
        all(
            c("sample", "locus", "allele", "allele_no", "reads") %in% names(gen)
        )
    )
    # error if malformed variant names
    var_mal <- rename(var_red, whatever = variant)
    expect_error(
        genotype(var_mal)
    )
})

# ploidy levels returns expected genotypes
test_that("'ploidy' parameter returns expected genotypes", {
    # gen with diff ploidy levels
    gen_hap <- genotype(var_red, ploidy = 1)
    gen_dip <- genotype(var_red, ploidy = 2)
    gen_pol <- genotype(var_red, ploidy = "poly")
    # other values returns error
    expect_error(
        genotype(var_red, ploidy = "whatever")
    )
    expect_error(
        genotype(var_red, ploidy = 3)
    )
    # haploid gen has one allele per locus
    expect_true(
        all(gen_hap$allele_no == 1)
    )
    # diploid gen has two allele per locus
    expect_true(
        max(gen_dip$allele_no) == 2
    )
})

# ADt sets hemizosity threshold as expected
test_that("ADt sets threshold for hemizygous correctly", {
    # select locus x sample for which there is only one variant, thus,
    # impliying homozygote or hemizygote
    var_h <- var_red[!duplicated(var_red$sample), ]
    gen_min <- genotype(var_h,
        ADt = min(var_h$reads) - 1
    )
    gen_more <- genotype(var_h,
        ADt = min(var_h$reads) * 2
    )

    # expect all variants are considered to be in homozygosis when.
    #   ADt is below number or reads.
    expect_equal(nrow(gen_min) / 2, nrow(var_h))
    expect_true(nrow(gen_min) > nrow(gen_more))
})
