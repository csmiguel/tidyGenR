test_that("BioString object is generated correctly", {
    data("genotypes")
    # the filtering from the fasta header constructor works ok
    seqs1 <-
        tidy2sequences(
            td = genotypes,
            fasta_header = "{locus}_{allele}",
            filename = FALSE
        )
    ndiff <- select(genotypes, locus, allele) |>
        distinct() |>
        nrow()
    expect_equal(length(seqs1), ndiff)
    # Biostrings is produced
    expect_equal(as.character(class(seqs1)), "DNAStringSet")
    # sequences are named according to fasta header
    ndiff_names <-
        select(genotypes, locus, allele) |>
        distinct() |>
        unite(sname, locus, allele) |>
        pull(sname)
    expect_equal(names(seqs1), ndiff_names)
})

test_that("FASTA is produced correctly:", {
    data("genotypes")
    # write FASTA to file
  p_seqs <- file.path(tempdir(), "test.fasta")
    seqs1 <-
        tidy2sequences(
            td = genotypes,
            fasta_header = "{locus}_{allele}",
            filename = p_seqs
        )

    seqs2 <- Biostrings::readDNAStringSet(p_seqs)
    # read FASTA is similar to BS output
    expect_equal(seqs1, seqs2)
})
