#' Evaluate alphanumeric strings
#'
#' Evaluates if all characters in a string are alphanumeric
#'
#' @param x String.
#' @returns T/F if all characters in a string are/not alphanumeric.
#' @noRd
#' @examples
#' x <- c("a", "bbb")
#' y <- c("a")
#' z <- c("a", "b_")
#' is.aphanumeric(x) # T
#' is.aphanumeric(y) # T
#' is.aphanumeric(z) # F
is_alphanumeric <- function(x) {
    if (!is.character(x)) {
        stop("x has to be of class character")
    }
    h <- strsplit(x = paste0(x, collapse = ""), split = "")
    hh <- vapply(h, function(x) all(grepl("[a-zA-Z0-9]", x)), logical(1))
    all(hh)
}

#' Check if any element in a character vector is not a DNA sequence
#'
#' If a character string read by DNAString is non a DNA sequence it will
#' return an error.
#' @return TRUE if any element is cannot be converted to DNAString;
#' FALSE if all elements can be converted to DNAString.
#' @param vector_seqs character vector with DNA sequences.
#' @noRd
any_non_DNAString <- function(vector_seqs) {
    any(
        vapply(vector_seqs, function(z) {
            is(try(DNAString(z), silent = TRUE), "try-error")
        }, logical(1))
    )
}

#' Check data frame with primers
#'
#' Evaluates the format of the data frame containing the name of the locus and
#' the primer sequences used. All names and sequences must be alphanumeric
#' strings, and three variables present: 1. locus, alphanumeric only loci
#' names, 2. fw, 5'->3' forward primer sequences, 3. rv, 5'->3'
#' reverse primer sequence.
#'
#' @param x  dataframe.
#' @returns If format is not correct it stops execution and returns
#' related error.
#' @noRd
check_primers <- function(x) {
    if (!"data.frame" %in% class(x) | !all(names(x) == c("locus", "fw", "rv"))) {
        stop(
            "'x' has to be a dataframe with 3 columns: 1. 'locus', ",
            "alphanumeric only loci names, 2. 'fw', 5'->3' forward ",
            "primer sequences, 3. 'rv', 5'->3' reverse primer sequence."
        )
    }
    if (!all(vapply(x$locus, is_alphanumeric, logical(1)))) {
        stop(
            "field 'locus' in dataframe has one or multiple strings ",
            "containing non-alphanumeric characters."
        )
    }
    if (any_non_DNAString(x$fw)) {
        stop(
            "field 'fw' in dataframe has one or multiple strings not ",
            "convertable to DNAString."
        )
    }
    if (any_non_DNAString(x$rv)) {
        stop(
            "field 'rv' in dataframe has one or multiple strings not ",
            "convertable to DNAString."
        )
    }
    message(
        "Dataframe with primers has expected columns with ",
        "loci names and primer sequences."
    )
}

#' Evaluate if all demultiplexed-by-locus-and-sample sequencing files
#'  meet the required format for
#' input in  the function 'trunc_amp'. That is: "sample.locus.\[1|2\].fastq.gz".
#' Evaluates if all characters in a string are alphanumeric
#'
#' @param ... other arguments from the environment.
#' @returns A list of \[1\] samples and \[2\] loci detected.
#' Stops if sample names do not conform the required format.
#' @noRd
#' @rdname truncate
check_names_demultiplexed <- function(in_dir,
                                      fw_pattern,
                                      rv_pattern = NULL) {
    fw_files <- sort(list.files(in_dir,
        pattern = fw_pattern,
        full.names = TRUE
    ))
    rv_files <- sort(list.files(in_dir,
        pattern = rv_pattern,
        full.names = TRUE
    ))
    # check if names meet the format "sample.locus.[1|2].fastq.gz".
    # list with names for forward reads
    lnamesf <-
        strsplit(basename(fw_files), "\\.")
    lnames <- lnamesf
    # list with names for reverse reads, if present
    if (length(rv_files) > 0) {
        lnamesr <-
            strsplit(basename(rv_files), "\\.")
        lnames <- c(lnamesf, lnamesr)
    }

    # print samples names detected
    snames <-
        sort(
            unique(
                vapply(lnames, function(x) x[[1]], character(1))
            )
        )
    snames_str <- paste(snames, collapse = ", ")
    message("\nSamples detected:\n ", snames_str)
    # print loci detected
    locinames <-
        sort(
            unique(
                vapply(lnames, function(x) x[[2]], character(1))
            )
        )
    locinames_str <- paste(locinames, collapse = ", ")
    message("\nLoci detected: \n", locinames_str)
    # test of all files read have alphanumeric sample name and locus
    # name in their names.
    nameformat_t <-
        vapply(lnames, function(x) {
            is_alphanumeric(x[1]) & is_alphanumeric(x[2])
        }, logical(1))

    if (!all(nameformat_t)) {
        stop(
            "Sequencing files matching patterns do not seem to ",
            "conform the required naming format 'sample.locus.[1|2].fastq.gz'"
        )
    } else if (all(nameformat_t)) {
        message(
            "\nSequencing files matching patterns seem to conform ",
            "the required naming format 'sample.locus.[1|2].fastq.gz'\n"
        )
    }

    return(list(loci = locinames, samples = snames))
}

#' Check F/R files have a match
#'
#' Checks that for each forward file there is a file with the same
#' name with R reads.
#'
#' @param fw_files character vector with path to forward files
#' @param rv_files character vector with path to reverse files
#' @param sample_locus Patterns to extract from FASTQ file names.
#' Group 1 captures
#' sample name and group 2 captures locus name.
#' (DEFAULT: `(^[a-zA-Z0-9]*)_([a-zA-Z0-9]*)`).
#' `^[a-zA-Z0-9]*_[a-zA-Z0-9]*` will extract 'sample_locus'
#' @param allpaired message to return if all files match.
#' @param not_paired "stop" message to return if F/R files do not match.
#' @returns message or stops R process. Matching substring in F/R reads
#'  are printed.
#' @noRd
#' @examples
#' check_fr_files(
#'     fw_files = filtFs,
#'     rv_files = filtRs,
#'     sample_locus = "^[a-zA-Z0-9]*_[a-zA-Z0-9]*",
#'     allpaired = "All F/R files match",
#'     not_paired = "Not all F/R files match"
#' )
check_fr_files <- function(fw_files,
                           rv_files,
                           sample_locus,
                           allpaired = "All F/R file names match.",
                           not_paired = "Not all F/R file names match.") {
  if(identical(fw_files, rv_files)) {
    stop("F and R file names provided are the same")
  }
    h <- basename(fw_files)
    hh <- sort(str_extract(h, sample_locus))
    hr <- basename(rv_files)
    hhr <- sort(str_extract(hr, sample_locus))
    if (length(hh) == 0 || length(hhr) == 0) {
        message("F or R have not been detected.")
    } else if (!all(hh == hhr)) {
        stop(not_paired)
    } else {
        message(allpaired)

        hh_files <- paste(hh, collapse = ", ")
        message("Files with F/R reads:\n", hh_files, "\n")
    }
}

#' Guess ploidy
#'
#' Tries to guess ploidy from tidy genotypes based on the maximum
#'  number of alleles observed for
#' any sample~locus combination.
#'
#' @param gen tidy genotypes.
#' @returns integer. Ploidy.
#' @noRd
guess_ploidy <- function(gen) {
    val <- plyr::daply(
        gen, ~ locus + sample,
        nrow
    ) |>
        as.numeric()
    val <- val[!is.na(val)]
    if (length(val) > 1) {
        warning("Guessed ploidy is likely ", max(val))
    }
    return(max(val))
}
