#' Truncate reads
#'
#' Truncates reads using 'filterAndTrim()' from DADA2. It automatically detects
#' single-end or paired-end sequencing files.
#'
#' @param loci Character vector with loci. If NULL, all loci detected by
#' 'check_names_demultiplexed()' are processed.
#' @param in_dir Path to folder with demultiplexed reads "sample.locus.
#' \[1|2\].fastq.gz".
#' @param fw_pattern Pattern matching specific extension of forward reads.
#' @param rv_pattern Pattern matching specific extension of reverse reads.
#' @param trunc_fr Flexible argument that sets truncation length for F
#' (and R) reads.
#'  - If 'numeric' of length == 1, trunFR is used for F and optionally R
#' reads, trunc_f = trunc_r.
#'  - If 'numeric' of length == 2, the first and second positions are used
#' as truncation lengths in F and R reads, accross all loci.
#'  - If it is a 'dataframe' with 'colnames == c("locus", "trunc_f", "trunc_r")
#' ', locus-spefic truncation lengths are applied from 'trunc_f' and 'trunc_r'
#'  columns. If 'se', 'trunc_r' values are ignored.
#' @param outdir Path where filtered reads are written. Created if it does
#'  not exist.
#' @param filt_name Pattern to append to FASTA names: '{sample}_{locus}{filt_name}'.
#' @param mode_trun 'se': single-end; 'pe', paired-end.
#' @param multithread T/F, see 'filterAndTrim()'.
#' @param max_ee Maximum expected errors. See 'filterAndTrim()'.
#' @param trunc_q Trim low-quality bases from 3' end.
#'  See 'filterAndTrim()'.
#' @return List with dataframes of number of reads before and after truncation.
#'  It writes truncated sequence to "outdir"
#' @name truncate
#' @examples
#' dem <-
#'  system.file("extdata", "demultiplexed", package = "tidyGenR")
#' trunc_amp(
#'     mode_trun = "pe",
#'     in_dir = dem,
#'     fw_pattern = "1.fastq.gz",
#'     rv_pattern = "2.fastq.gz",
#'     trunc_fr = c(250, 180),
#'     max_ee = c(3, 3)
#' )
#' # example 2:
#' data("trunc_fr")
#' trunc_amp(
#'     mode_trun = "se",
#'     loci = c("chrna9", "nfkbia"),
#'     in_dir = dem,
#'     fw_pattern = "1.fastq.gz",
#'     rv_pattern = "2.fastq.gz",
#'     trunc_fr = trunc_fr,
#'     max_ee = 3,
#'     trunc_q = 2
#' )
#' @export
trunc_amp <- function(loci = NULL,
                      in_dir, trunc_fr,
                      fw_pattern, rv_pattern = NULL,
                      outdir = "truncated",
                      filt_name = "_F_filt.fastq.gz",
                      mode_trun = "pe", multithread = FALSE,
                      max_ee = 3, trunc_q = 2) {
    stopifnot(dir.exists(in_dir))
    z <-
        check_names_demultiplexed(
            in_dir = in_dir,
            fw_pattern = fw_pattern,
            rv_pattern = rv_pattern
        )
    if (is.null(loci)) {
        loci <- z$loci
    }
    filtering_output <-
        lapply(loci, function(locus) {
            message("\nTruncating ", locus, "\n")
            # forward fastq files for a given locus
            fws <-
                get_in_and_out_filenames(
                    locus = locus,
                    pattern = fw_pattern,
                    in_dir = in_dir,
                    outdir = outdir,
                    filt_pattern = filt_name
                )
            # check there are input reads
            if (length(fws$infiles) == 0)
              stop("No files have been found for ", locus, " locus")
            # truncation lengths
            trl <-
                set_tunc_fr(
                    trunc_fr = trunc_fr,
                    locus = locus,
                    mode_trun = mode_trun
                )
            # check well-formed truncation length
            if(!is.numeric(trl))
              stop("Truncation length/s for ", locus, " is not numeric.")
            if (mode_trun == "se") {
                out <-
                    trunc_se(
                        fws = fws, locus = locus,
                        trl = trl, max_ee = max_ee,
                        trunc_q = trunc_q,
                        multithread = multithread
                    )
            } else if (mode_trun == "pe") {
                out <-
                    trunc_pe(
                        fws = fws, locus = locus, rv_pattern = rv_pattern,
                        in_dir = in_dir, outdir = outdir,
                        trl = trl, max_ee = max_ee,
                        trunc_q = trunc_q,
                        multithread = multithread
                    )
            }

            if (sum(out[, 2]) == 0) {
                message("\nFiltering removed all reads for", locus, ".\n")
            }
            rownames(out) <- fws$snames
            return(out)
        })
    setNames(filtering_output, loci)
}

#' helper functions to trunc_amp: get truncation lengths
#' @param locus Specific locus.
#' @rdname truncate
#' @noRd
set_tunc_fr <- function(trunc_fr, locus, mode_trun) {
    trunc_fr_locusi <- numeric(2L)
    if ("data.frame" %in% class(trunc_fr)) {
        stopifnot(identical(names(trunc_fr), c("locus", "trunc_f", "trunc_r")))
        stopifnot(locus %in% trunc_fr$locus)
        trunc_fr_locusi[1] <- trunc_fr$trunc_f[trunc_fr$locus == locus]
        trunc_fr_locusi[2] <- trunc_fr$trunc_r[trunc_fr$locus == locus]
    } else if (is.numeric(trunc_fr)) {
        if (length(trunc_fr) == 1 && mode_trun == "se") {
            trunc_fr_locusi[1] <- trunc_fr
        }
        if (length(trunc_fr) == 1 && mode_trun == "pe") {
            trunc_fr_locusi[1] <- trunc_fr
            trunc_fr_locusi[2] <- trunc_fr
        }
        if (length(trunc_fr) == 2) {
            trunc_fr_locusi[1] <- trunc_fr[1]
            trunc_fr_locusi[2] <- trunc_fr[2]
        }
    } else {
        stop("'trunc_fr' is not valid.")
    }
    return(trunc_fr_locusi)
}

#' helper functions to trunc_amp: get truncation lengths
#' @param filt_pattern Appended suffix to filtered FASTQ. It has been
#' fixed to `_[F|R]_filt.fastq.gz`.
#' @rdname truncate
#' @returns A list with 3 elements:
#'  - *infiles*: character vector with full path to input FASTQ.
#'  - *outfiles*: character vector with full path to filtered FASTQ.
#'  - *snames*: character vector to sample names.
#' @noRd
get_in_and_out_filenames <- function(locus,
                                     pattern,
                                     in_dir,
                                     outdir,
                                     filt_pattern) {
    p <- paste0(".", locus, ".", pattern)
    # list files
    f <- sort(list.files(in_dir,
        pattern = p,
        full.names = TRUE
    ))
    # samples from the fastq files
    sample_names <-
        vapply(strsplit(basename(f), "\\."), function(x) x[1], character(1))
    # path to filtered
    filt <-
        file.path(
            outdir,
            paste0(
                sample_names, "_",
                locus,
                filt_pattern
            )
        )
    return(list(infiles = f, outfiles = filt, snames = sample_names))
}

#' helper functions to trunc_amp:
#' truncate using single-end method for a given LOCUS.
#' @noRd
trunc_pe <- function(fws, locus, in_dir, rv_pattern, trl,
                     outdir, multithread, max_ee, trunc_q) {
    stopifnot(length(fws$infiles) > 0)
    rvs <-
        get_in_and_out_filenames(
            locus = locus,
            pattern = rv_pattern,
            in_dir = in_dir,
            outdir = outdir,
            filt_pattern = "_R_filt.fastq.gz"
        )
    stopifnot(length(rvs$infiles) == length(fws$infiles))
    out <- filterAndTrim(fws$infiles, fws$outfiles,
        rvs$infiles, rvs$outfiles,
        maxN = 0,
        compress = TRUE,
        truncLen = trl,
        multithread = multithread,
        maxEE = max_ee,
        truncQ = trunc_q
    )
    message(paste(
        "\nForwards reads of", locus,
        "truncated at", trl[1], "nt."
    ))
    message(paste(
        "\nReverse reads of", locus,
        "truncated at", trl[2], "nt.\n"
    ))
    return(out)
}

#' helper functions to trunc_amp:
#' truncate using paired-end method for a given LOCUS.
#' @noRd
trunc_se <- function(fws, locus, trl,
                     multithread, max_ee, trunc_q) {
    stopifnot(length(max_ee) == 1)
    out <- filterAndTrim(fws$infiles, fws$outfiles,
        maxN = 0,
        compress = TRUE,
        truncLen = trl[1],
        multithread = multithread,
        maxEE = max_ee,
        truncQ = trunc_q
    )
    message(paste(
        "\nForwards reads of", locus, "truncated at",
        trl[1], "nt.\n"
    ))
    return(out)
}
