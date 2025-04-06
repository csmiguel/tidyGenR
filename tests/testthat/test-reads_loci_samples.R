truncated <-
  system.file("extdata", "truncated", package = "tidyGenR")

# 'path' accepts a folder location or
test_that("Expected default dataframes are produced when using a path to folder", {
    df <-
        reads_loci_samples(
            path = truncated,
            pattern_fq = "F_filt.fastq.gz"
        )
    df2 <-
        reads_loci_samples(
            path = truncated,
            pattern_fq = "F_filt.fastq.gz",
            all.variants = "whatever",
            var_id = "whatever"
        )
    expect_equal(dim(df), c(3, 4))
    # all.variants and var_id are ignored when setting path to folder
    expect_identical(df, df2)
})

test_that("the expected dataframe are output when path is a data.frame", {
    data("variants")

    df_default <-
        reads_loci_samples(path = variants)

    df_var_true <-
        reads_loci_samples(
            path = variants,
            all.variants = TRUE
        )
    df_var_var <-
        reads_loci_samples(
            path = variants,
            all.variants = TRUE,
            var_id = "variant"
        )

    # same dimensions when change var_id
    expect_equal(dim(df_var_true), dim(df_var_var))
    # greater number of cols when all.variants is TRUE
    expect_true(ncol(df_default) < ncol(df_var_true))
    # default formation from colnames matches var_id
    in_var_names <- variants$md5
    out_var_names <- stringr::str_remove(colnames(df_var_true), "^.*_")
    expect_true(all(in_var_names %in% out_var_names))
    # custom formation from colnames eg, variants, matches var_id
    in_var_names <- variants$variant
    out_var_names <- stringr::str_remove(colnames(df_var_var), "^.*_")
    expect_true(all(in_var_names %in% out_var_names))
})
