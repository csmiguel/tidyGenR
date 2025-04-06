test_that("col and row filtering in msa returns expected nseqs and width", {
    data("msa")
    nseqs <- 2
    ncol_msa <- 50
    msa_filt <- filt_msa(msa, 1:nseqs, c(1:ncol_msa))
    expect_equal(length(msa_filt), nseqs)
    expect_equal(unique(Biostrings::width(msa_filt)), ncol_msa)
})
