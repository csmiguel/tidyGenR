truncated <-
  system.file("extdata", "truncated", package = "tidyGenR")

test_that("Default parameters work for single end:", {
    variants <- variant_call(in_dir = truncated)
    expect_equal(colnames(variants), c(
        "sample", "locus", "variant", "reads",
        "nt", "md5", "sequence"
    ))
})

test_that("loess_err_mod4 returns no error", {
    variants <-
        variant_call(
            in_dir = truncated,
            error_function = loess_err_mod4
        )
    expect_equal(colnames(variants), c(
        "sample", "locus", "variant", "reads",
        "nt", "md5", "sequence"
    ))
})

test_that(paste(
    "setting a known locus name calls variants for",
    "that given locus and return expected sample and locus names"
), {
    locus <- "chrna9"
    variants <-
        variant_call(
            loci = locus,
            in_dir = truncated
        )
    expect_equal(nrow(variants), 5)
    expect_equal(as.character(unique(variants$locus)), locus)
    expect_equal(levels(as.factor(variants$sample)),
                 c("BOR1061", "BOR1063", "BOR1069"))
})

test_that("single end and paired-end mode work", {
    locus <- "rogdi"
    single_end <-
        variant_call(
            loci = locus,
            in_dir = truncated
        )
    paired_end <-
        variant_call(
            loci = locus,
            in_dir = truncated,
            rv_pattern = "_R_filt.fastq.gz"
        )
    expect_equal(unique(single_end$nt), 270L)
    expect_equal(unique(paired_end$nt), 398L)
})

test_that("'c_unmerge' does not have an effect in single-end mode or when F/R reads overlap", {
    locus <- "chrna9" # all reads overlap
    # paired end
    paired_end_merged <-
        variant_call(
            loci = locus,
            in_dir = truncated,
            rv_pattern = "_R_filt.fastq.gz",
            c_unmerged = TRUE
        )
    paired_end_unmerged <-
        variant_call(
            loci = locus,
            in_dir = truncated,
            rv_pattern = "_R_filt.fastq.gz",
            c_unmerged = FALSE
        )
    expect_equal(paired_end_merged, paired_end_unmerged)
    # single end
    single_end_merged <-
        variant_call(
            loci = locus,
            in_dir = truncated,
            rv_pattern = NULL,
            c_unmerged = TRUE
        )
    single_end_unmerged <-
        variant_call(
            loci = locus,
            in_dir = truncated,
            rv_pattern = NULL,
            c_unmerged = FALSE
        )
    expect_equal(single_end_merged, single_end_unmerged)
})

test_that("samples are pooled when pool == T, and non pooled when pool == F", {
    locus <- "chrna9"
    # pooled
    out_pooled <- capture_output(variant_call(
        loci = locus,
        in_dir = truncated,
        rv_pattern = "_R_filt.fastq.gz",
        pool = TRUE
    ))
    expect_true(grepl("samples were pooled", out_pooled))
    # non pooled
    out_non_pooled <- capture_output(variant_call(
        loci = locus,
        in_dir = truncated,
        rv_pattern = "_R_filt.fastq.gz",
        pool = FALSE
    ))
    expect_false(grepl("samples were pooled", out_non_pooled))
})


# test_that("multithread is activate when multithread == T", {
# })

test_that("chimera removal can be tunned with different methods or not run", {
    locus <- "chrna9"
    # chim_rm == FALSE
    chim_rm_false <-
        capture_messages(
            variant_call(
                loci = locus,
                in_dir = truncated,
                rv_pattern = "_R_filt.fastq.gz",
                chim_rm = FALSE
            )
        )
    expect_equal(
        tail(chim_rm_false, 1),
        "Chimera removal step has been skipped.\n"
    )
    # chim_rm == consensus
    chim_rm_consensus <-
        capture_messages(
            variant_call(
                loci = locus,
                in_dir = truncated,
                rv_pattern = "_R_filt.fastq.gz",
                chim_rm = "consensus"
            )
        )
    expect_true(grepl(pattern = "bimera", tail(chim_rm_consensus, 1)))
    # chim_rm assumes other arguments as in removeBimeraDenovo()
    chim_rm_pooled <-
        capture_messages(
            variant_call(
                loci = locus,
                in_dir = truncated,
                rv_pattern = "_R_filt.fastq.gz",
                chim_rm = "pooled"
            )
        )
    expect_true(grepl(pattern = "bimera", tail(chim_rm_pooled, 1)))
})

test_that("AD and MAF are working", {
    locus <- "chrna9"
    ad <- 25
    var_unfiltered <-
        variant_call(
            loci = locus,
            in_dir = truncated
        )
    var_filtered <-
        variant_call(
            loci = locus, maf = 0.5, ad = ad,
            in_dir = truncated
        )
    expect_equal(nrow(var_unfiltered), 5)
    expect_equal(nrow(var_filtered), 3)
    expect_gt(max(var_filtered$reads), ad)
})



# unfinished
test_that("tunning  omega a has an effect on the number of variants determined.", {
    # sensitivity to omega_r
    variants_omegar_high <-
        variant_call(
            loci = "chrna9",
            rv_pattern = "_R_filt.fastq.gz",
            in_dir = truncated,
            omega_a_r = 0.1
        )
    variants_omegar_low <-
        variant_call(
            loci = "chrna9",
            rv_pattern = "_R_filt.fastq.gz",
            in_dir = truncated
        )
    expect_true(nrow(variants_omegar_high) > nrow(variants_omegar_low))
})

