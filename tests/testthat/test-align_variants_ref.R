data("variants")
ref_al <-
  system.file("extdata/reference_alleles.fasta", package = "tidyGenR")

test_that("multiplication works", {
  al <- align_variants_ref(variants, ref_al)
  al_filt <- al[al$Score > 100, ]
  testthat::expect_true(nrow(al_filt) > 10)
})
