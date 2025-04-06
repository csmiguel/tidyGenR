data("genotypes")

test_that("Hemizygous are removed for diploid data", {
    gen2 <- dplyr::slice(genotypes, 1:5)
    gen3 <- suppressWarnings(remove_hemizygotes(gen2))
    expect_true(nrow(gen3) < nrow(gen2))
})

test_that("returns error for haploid data", {
    gen2 <- dplyr::slice(genotypes, c(1, 3, 5))
    expect_error(suppressWarnings(remove_hemizygotes(gen2)),
        regexp = "Genotypes have to be diploid."
    )
})
