test_that("multiplication works", {
    data("variants")
    seqs <- variants$sequence[1:5]
    bs <- view_msa(seqs, browse = FALSE, out_bs = TRUE)
    expect_equal(as.character(class(bs)), "DNAStringSet")
})
