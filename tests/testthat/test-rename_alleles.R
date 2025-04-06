test_that("default replace_alleles works with input ref dataframe and fasta file", {
    data("genotypes")
    ref_al <-
        system.file("extdata/reference_alleles.fasta", package = "tidyGenR")
    # toy df with alleles named as integers
    df_toy_ref <-
        gen_tidy2integers(genotypes) |>
        select(allele, sequence) |>
        mutate(allele = as.character(allele)) |>
        distinct()
    # rename alleles with integers
    df_renamed_df <-
        rename_alleles(
            df_alleles = genotypes,
            ref_alleles = df_toy_ref
        )
    # way B to rename alleles by integers
    df_renamed_int <-
        gen_tidy2integers(genotypes) |>
        mutate(allele = as.character(allele)) |>
        arrange(sample, locus)

    # check renaming alleles with 'rename_alleles' returns
    # same data as 'gen_tidy2integers'
    expect_equal(
        head(data.frame(df_renamed_df)),
        head(data.frame(df_renamed_int))
    )

    # rename based on a FASTA file
    df_renamed_fa <-
        rename_alleles(
            df_alleles = genotypes,
            ref_alleles = ref_al
        )
    expect_equal(genotypes$sequence, df_renamed_fa$sequence)
    expect_true(all(genotypes$sequence == df_renamed_fa$sequence))
    expect_false(all(genotypes$allele == df_renamed_fa$allele))
})

# test 'replacement' works

test_that("'replacement' returns 2-column df", {
    data("genotypes")
    ref_al <-
        system.file("extdata/reference_alleles.fasta", package = "tidyGenR")
    rep_al <-
        rename_alleles(
            df_alleles = genotypes,
            ref_alleles = ref_al,
            replacements = TRUE
        )
    expect_equal(names(rep_al), c("locus", "old_name", "new_name"))
})

# error if mandatory fields are not present
test_that("error when 'allele' or 'sequence' fields are not present", {
    data("genotypes")
    ref_al <-
        system.file("extdata/reference_alleles.fasta", package = "tidyGenR")
    geno2 <- rename(genotypes, sekuence = sequence)
    expect_error(rename_alleles(
        df_alleles = geno2,
        ref_alleles = ref_al
    ))
})
