#' Variant calling using DADA2 ASV approach
#'
#' @description
#' - 'variant_call_dada()' Call variants for a single locus using DADA2.
#' It is nested in 'variant_call'.
#' @details
#'  - 'variant_call_dada()' Allows single-end and paired-end data.
#' 'c_unmerged' will trigger the 'justConcatenate' argument in 'mergePairs',
#' and 10 N's will be used to concatenate non-overlapping F and R reads.
#' Use 'concateNotMerged' carefully, as it will generate artificial variant
#'  sequences. Default is deactivated.
#'
#' @param locus Locus name.
#' @param in_dir Path to folder with truncated files.
#' @param fw_pattern Pattern matching files with F reads.
#' @param rv_pattern Pattern matching files with R reads. If left NULL,
#' single-end sequencing will be assumed.
#' @param sample_locus Patterns to extract from FASTQ file names.
#' Group 1 captures
#' sample name and group 2 captures locus name.
#' (DEFAULT: `(^[a-zA-Z0-9]*)_([a-zA-Z0-9]*)`).
#' `^[a-zA-Z0-9]*_[a-zA-Z0-9]*` will extract 'sample_locus'
#' from default naming convention `sample_locus_[F|R]_fastq.gz`.
#' @param c_unmerged F/R sequences that were not merged in mergePairs
#' are concatenated using a stretch of 10 N's.
#' @param pool Passed to 'dada()'. Denoising is done in pooled samples
#' (T) or by sample (F).
#' @param error_function Use default 'loessErrfun()' for regular
#' Illumina quality codification and 'loess_err_mod4()' for binned NovaSeq
#' qualities. Passed to 'dada(errorEstimationFunction)'.
#' @param multithread T/F, passed to 'dada(multithread)' and
#' 'learnErrors(multithread)'.
#' @param chim_rm If FALSE, no chimera removal is performed. If chim_rm  is
#' "character", it is passed to 'removeBimeraDenovo(method)'.
#' @param omega_a_f "OMEGA_A" passed to 'dada' in forward reads
#'  (Default: 'getDadaOpt()$OMEGA_A').
#' @param omega_a_r "OMEGA_A" passed to 'dada' in reverse reads
#' (Default: 'getDadaOpt()$OMEGA_A').
#' @param band_size "BAND_SIZE" passed to 'dada'
#' (DEFAULT: 'getDadaOpt()$BAND_SIZE').
#' @return
#' - 'variant_call_dada()' Numeric array. Column names are variant sequences;
#' row names are sample names and each number in the array indicate
#' the number of reads supporting variant 'i' in sample 'j'.
#' @name variant_calling_dada
#' @rdname variant_calling_dada
variant_call_dada <- function(
    locus,
    in_dir,
    fw_pattern = "_F_filt.fastq.gz",
    rv_pattern = NULL,
    sample_locus = "(^[a-zA-Z0-9]*)_([a-zA-Z0-9]*)",
    c_unmerged = FALSE,
    pool = FALSE,
    error_function = loessErrfun,
    multithread = FALSE,
    chim_rm = "consensus",
    omega_a_f = getDadaOpt()$OMEGA_A,
    omega_a_r = getDadaOpt()$OMEGA_A,
    band_size = getDadaOpt()$BAND_SIZE) {
    # 0. get forward and (reverse files)
    message("\nCalling variants for ", locus, "\n")
    # list of fastq files for a given locus
    fw_filt_all <-
        sort(list.files(in_dir,
            pattern = fw_pattern,
            full.names = TRUE
        ))
    # filter by locus
    lnames <-
        str_extract(basename(fw_filt_all), sample_locus, group = 2)
    fw_fq <- sort(fw_filt_all[which(lnames == locus)])
    # sample names
    snames <-
        str_extract(basename(fw_fq), sample_locus, group = 1)
    # For se and pe:
    err_fw <-
        learnErrors(fw_fq,
            multithread = multithread,
            errorEstimationFunction = error_function
        )
    dada_fw <-
        dada(fw_fq,
            err = err_fw,
            multithread = multithread,
            pool = pool,
            verbose = TRUE,
            OMEGA_A = omega_a_f,
            BAND_SIZE = band_size
        ) |>
      dada2list(names = basename(fw_fq))
    seqtab <-
        makeSequenceTable(dada_fw)
    # check at least a sample has reads after dada
    z <-
        vapply(dada_fw, function(x) sum(x$clustering$abundance), numeric(1))
    if (!sum(z) > 0) {
        stop("No sequences after dada.")
    }
    if (!is.null(rv_pattern)) {
        rv_filt_all <-
            sort(list.files(in_dir,
                pattern = rv_pattern,
                full.names = TRUE
            ))
        stopifnot(length(rv_filt_all) > 0)
        # filter by locus
        lnames <-
            str_extract(basename(rv_filt_all), sample_locus, group = 2)
        rv_fq <- sort(rv_filt_all[which(lnames == locus)])
        # check there are paired F/R files
        check_fr_files(
            fw_files = fw_fq,
            rv_files = rv_fq,
            sample_locus = sample_locus
        )
        err_rv <-
            learnErrors(rv_fq,
                multithread = multithread,
                errorEstimationFunction = error_function
            )
        dada_rv <-
            dada(rv_fq,
                err = err_rv,
                multithread = multithread,
                pool = pool,
                verbose = TRUE,
                OMEGA_A = omega_a_r,
                BAND_SIZE = band_size
            )
        dada_rv <- dada2list(dada_rv, names = basename(rv_fq))
        message("\nF and R reads found. Merging reads:\n")
        dadamerged <-
            mergePairs(dada_fw, fw_fq,
                dada_rv, rv_fq,
                verbose = TRUE, minOverlap = 10,
                maxMismatch = 0, trimOverhang = TRUE
            )
        dadamerged_list <- dada2list(dadamerged, names = basename(fw_fq))
        seqtab <-
            makeSequenceTable(dadamerged_list)
        samples_unmerged <-
            which(
                vapply(dadamerged_list, function(x) sum(x$abundance), numeric(1)) == 0
            )
        if (length(samples_unmerged) == 0) {
            message(
                "All samples have overlapping F/R reads that have been merged ",
                "successfully."
            )
        } else if (length(samples_unmerged) > 0) {
            message(
                "Sample/s ",
                paste(snames[samples_unmerged], collapse = " "),
                " reads has/have not been merged successfully."
            )
        }
        ### WARNING # reads might no merge if there are no overlapping regions.
        # There is an option to concatenate them 'c_unmerged = T'. It will
        # concatenate the reads for those samples for which the merging yielded 0
        # counts.
        if (c_unmerged && length(samples_unmerged) > 0) {
            # if any of the samples have reads wh
            message(
                "'c_unmerged' is activated and merging did not work for some",
                " samples. Concatenating F/R reads for ",
                paste(snames[samples_unmerged],
                    collapse = " "
                )
            )
            dadaconcat <-
                mergePairs(dada_fw[samples_unmerged], fw_fq[samples_unmerged],
                    dada_rv[samples_unmerged], rv_fq[samples_unmerged],
                    verbose = TRUE, minOverlap = 10,
                    maxMismatch = 0, trimOverhang = TRUE,
                    justConcatenate = TRUE
                )
            # merge datasets from merged reads and concatenated reads.
            # if only one sample is returned by mergePairs,
            # by default it will be data.frame and not a list.
            dadaconcat_list <- dada2list(dadaconcat, names(dada_fw)[samples_unmerged])

            dada_merged_concat <-
                c(
                    dadamerged_list[-samples_unmerged],
                    dadaconcat_list
                )
            seqtab <-
                makeSequenceTable(dada_merged_concat)
        }
    }
    rownames(seqtab) <-
        str_extract(rownames(seqtab), sample_locus, 1)
    # remove bimeras
    seqtab_no_chim <- chimera_removal(seqtab, chim_rm)
    return(seqtab_no_chim)
}


#' Loess fit to estimate error rates from transition counts in NovaSeq
#' binned qualities
#' @description
#' - 'loess_err_mod4()' Replaces default 'loessErrfun' when estimating error
#'  matrices for NovaSeq binned qualities.
#' @details
#' - 'loess_err_mod4()' This is a beta function which has been shown to work
#' in my experience and other users' experience.
#' See discussion in https://github.com/benjjneb/dada2/issues/1307
#' @param trans See loessErrfun for details.
#' @rdname variant_calling_dada
loess_err_mod4 <- function(trans) {
    qq <- as.numeric(colnames(trans))
    est <- matrix(0, nrow = 0, ncol = length(qq))
    for (nti in c("A", "C", "G", "T")) {
        for (ntj in c("A", "C", "G", "T")) {
            if (nti != ntj) {
                errs <- trans[paste0(nti, "2", ntj), ]
                tot <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
                rlogp <- log10((errs + 1) / tot)
                rlogp[is.infinite(rlogp)] <- NA
                df <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
                mod.lo <- loess(rlogp ~ q, df,
                    weights = log10(tot),
                    degree = 1, span = 0.95
                )
                pred <- predict(mod.lo, qq)
                maxrli <- max(which(!is.na(pred)))
                minrli <- min(which(!is.na(pred)))
                pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
                pred[seq_along(pred) < minrli] <- pred[[minrli]]
                est <- rbind(est, 10^pred)
            }
        }
    }

    # HACKY
    MAX_ERROR_RATE <- 0.25
    MIN_ERROR_RATE <- 1e-7
    est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
    est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE

    # enforce monotonicity
    # https://github.com/benjjneb/dada2/issues/791
    estorig <- est
    est <- est |>
        data.frame() |>
        mutate_all(list(case_when(
            . < X40 ~ X40,
            . >= X40 ~ .
        ))) |>
        as.matrix()
    rownames(est) <- rownames(estorig)
    colnames(est) <- colnames(estorig)

    # Expand the err matrix with the self-transition probs
    err <- rbind(
        1 - colSums(est[seq_len(3), ]), est[seq_len(3), ],
        est[4, ], 1 - colSums(est[4:6, ]), est[5:6, ],
        est[7:8, ], 1 - colSums(est[7:9, ]), est[9, ],
        est[10:12, ], 1 - colSums(est[10:12, ])
    )
    rownames(err) <- paste0(
        rep(c("A", "C", "G", "T"), each = 4),
        "2", c("A", "C", "G", "T")
    )
    colnames(err) <- colnames(trans)
    # Return
    return(err)
}

#' Remove chimeras
#' @description
#' - 'chimera_removal()' Wrapping function to remove chimeras
#' using 'removeBimeraDenovo()'.
#' @param seqtab Features table from dada2.
#' @return
#' - 'chimera_removal()' 'seqtab' like object (see dada2).
#' @rdname variant_calling_dada
chimera_removal <- function(seqtab, chim_rm) {
    if (dim(seqtab)[2] == 0) {
        no_chim <- NULL
    } else if (dim(seqtab)[2] > 0) {
        if (isFALSE(chim_rm)) {
            no_chim <- seqtab
            message("Chimera removal step has been skipped.")
        } else if (is.character(chim_rm)) {
            no_chim <-
                removeBimeraDenovo(
                    seqtab,
                    method = chim_rm,
                    verbose = TRUE
                )
        }
    }
    return(no_chim)
}

#' Convert dada2 output to list
#' @details
#' - 'dada2list()' DADA2 functions return a dataframe whenever
#' the length of the object
#' is == 1 or a list if length > 1. To smooth prevent from bugs, this function
#' standardizes the output of dada2 to lists when length == 1.
#' Unaltered when input is a list.
#' @param dada_object Object from [dada2::dada()], [dada2::mergePairs()], etc.
#' @param names Character vector to name the elements in the list. By default,
#'  DADA2 names are basenames from FAST(Q).
#' @return
#' - 'dada2list()' List with dada2 object/s.
#' @rdname variant_calling_dada
dada2list <- function(dada_object, names) {
    if (!"list" %in% class(dada_object)) {
        dada_object_l <- list(dada_object)
        names(dada_object_l) <- names
    } else {
        dada_object_l <- dada_object
    }
    return(dada_object_l)
}
