test_that("popart nex is well formed", {
    data("genotypes")
    x <-
        filter(genotypes, locus == "abcg8") %>%
        mutate(sname = paste0("seq_", seq_len(nrow(.))))
    msa <-
        tidy2sequences(x, fasta_header = "{sname}") %>%
        AlignSeqs()
    temp_file <- file.path(tempdir(), "test_popart.nex")
    out_popart(msa, x,
        outnex = temp_file,
        sname = "sname",
        xgroups = "sample",
        blocks = "TRAITS"
    )
    expect_snapshot_file(temp_file, "popart_test_nex")

    # adding DATA produces no error
    out_popart(msa, x,
        outnex = temp_file,
        sname = "sname",
        xgroups = "sample",
        blocks = c("TRAITS", "DATA")
    )
})
