test_that("multiplication works", {
    # create empty fastq.gz
    fq0 <-
        file.path(tempdir(), "empty_file.fastq")
    system(paste("touch", fq0))

    # create fastq with many seqs
    data("msa")
    fq_many <-
        file.path(tempdir(), "msa_file.fastq")
    Biostrings::writeXStringSet(msa, filepath = fq_many, format = "fasta")
    # remove empty fastq
    remove_poor_fastq(tempdir(),
        pattern = "fastq", min_reads = 1
    )
    # only empty file is removed
    expect_true(
        file.exists(fq_many)
    )
    expect_false(
        file.exists(fq0)
    )
})
