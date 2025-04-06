#' Variant calling
#'
#' @description
#' 'variant_call' Call variants for multiple loci.
#' @details
#' Allows single-end and paired-end data. Be careful with the use of
#' 'c_unmerged'.
#' It  will trigger the 'justConcatenate' argument in 'mergePairs', and
#' 10 N's will be used to concatenate non-overlapping F and R reads.
#' Use 'c_unmerged' carefully, as it will generate artificial
#' variant sequences.
#' Default is deactivated.
#' Variants are reformatted and to a "tibble" and filtered according to 'maf'
#' and 'ad', which are
#' added as attributes to output. Variables of output are 'sample', 'locus',
#' 'sequence' (DNA sequence), 'variant' (name of variant),
#' 'reads' (number of reads supporting that variant), 'nt' (sequence length),
#' 'md5' (md5 checksum).
#'
#' @param loci Character vector with loci to detect from 'sample_locus'
#' in sample names. If NULL, all loci detected according to the pattern in the
#' target directory are used.
#' @param in_dir Path to folder with truncated files.
#' @param fw_pattern Pattern matching files with F reads. It has to match
#' everything after locus name. eg see parenthesis in
#' "samplex_locus(_2_filt.fastq.gz)"
#' @param rv_pattern Pattern matching files with R reads. If left NULL,
#' single-end sequencing will be assumed.
#' @param sample_locus Patterns to extract from FASTQ file names.
#' Group 1 captures
#' sample name and group 2 captures locus name.
#' (DEFAULT: `(^[a-zA-Z0-9]*)_([a-zA-Z0-9]*)`).
#' `^[a-zA-Z0-9]*_[a-zA-Z0-9]*` will extract 'sample_locus'
#' from default naming convention `sample_locus_[F|R]_fastq.gz`.
#' @param c_unmerged F/R sequences that were not merged in mergePairs are
#'  concatenated using a stretch of 10 N's.
#' @param pool Passed to 'dada()'. Denoising is done in pooled samples
#' (T) or by sample (F).
#' @param error_function Use default 'loessErrfun' for regular Illumina quality
#'  codification and 'loess_err_mod4' for binned NovaSeq qualities.
#' @param multithread T/F, passed to 'multithread' in 'dada'
#' and 'learnErrors()'.
#' @param chim_rm If FALSE, no chimera removal is performed.
#' If == "character", it is passed to
#' 'method' in 'removeBimeraDenovo()' (DEFAULT = 'consensus').
#' @param maf Minimum Allele Frequency. Passed to filter_variants.
#' @param ad Allele Depth. Passed to filter_variants.
#' @param omega_a_f "OMEGA_A" passed to 'dada' in forward reads
#' (Default: getDadaOpt()$OMEGA_A).
#' @param omega_a_r "OMEGA_A" passed to 'dada' in reverse reads
#' (Default: getDadaOpt()$OMEGA_A).
#' @param band_size "BAND_SIZE" passed to 'dada'
#' (DEFAULT: getDadaOpt()$BAND_SIZE).
#' @return tidy tibble with locus, sample, sequence, variant (name),
#' nt(sequence length), md5.
#' @examples
#' # truncated fastq
#'truncated <-
#'  system.file("extdata", "truncated", package = "tidyGenR")
#' # variant calling
#' variant_call(in_dir = truncated)
#' @export
variant_call <- function(loci = NULL,
                         in_dir = NULL,
                         fw_pattern = "_F_filt.fastq.gz",
                         rv_pattern = NULL,
                         sample_locus = "(^[a-zA-Z0-9]*).([a-zA-Z0-9]*)",
                         c_unmerged = FALSE,
                         pool = FALSE,
                         error_function = loessErrfun,
                         multithread = FALSE,
                         chim_rm = "consensus",
                         ad = 1,
                         maf = 0,
                         omega_a_f = getDadaOpt()$OMEGA_A,
                         omega_a_r = getDadaOpt()$OMEGA_A,
                         band_size = getDadaOpt()$BAND_SIZE) {
    # block 0: get loci
    infil <-
        list.files(in_dir, pattern = fw_pattern)
    all_loc <-
        unique(str_extract(infil, pattern = sample_locus, group = 2))
    if (is.null(loci)) {
        loci <- all_loc
    } else if (is.character(loci)) {
        all_loc_true <- which(loci %in% all_loc)
        stopifnot(all(all_loc_true))
        loci <- loci[all_loc_true]
    } else {
        stop("'loci' is not from a expected 'class'.")
    }
    message(c(
        paste0(
            "A total of ",
            length(loci),
            " locus/loci have/has been detected in file(s) with the pattern '",
            fw_pattern,
            "' in the path '",
            in_dir, "':\n"
        ),
        paste(loci, collapse = " ")
    ))
    stopifnot(length(loci) > 0)
    if (is.null(rv_pattern)) {
        message("Trying to call variants in single end mode.")
    } else if (!is.null(rv_pattern)) {
        message("Trying to call variants in paired-end mode.")
    }
    # block 1: stream dada variant calling
    v1 <-
        lapply(loci, function(x) {
            variant_call_dada(
                locus = x, in_dir = in_dir,
                fw_pattern = fw_pattern, rv_pattern = rv_pattern,
                sample_locus = sample_locus, c_unmerged = c_unmerged,
                pool = pool, multithread = multithread, chim_rm = chim_rm,
                omega_a_f = omega_a_f,
                omega_a_r = omega_a_r, band_size = band_size
            )
        })
    names(v1) <- loci
    # remove loci with no variants
    # (they are returned as NULL elements in the list from variant_call_dada)
    v1 <- v1[!vapply(v1, is.null, logical(1))]
    # stop if no variant was found
    if (length(v1) == 0) {
        stop("No variant has been found for none of the loci nor samples.")
    }
    # second block: reformat variants, name alleles.
    # convert list to dataframe
    w <-
        ldply(v1, .id = "locus", function(x) {
            z <- rownames_to_column(data.frame(x), "sample")
            zz <- pivot_longer(z,
                cols = -sample,
                names_to = "sequence",
                values_to = "reads"
            )
        })
    ww <- filter_variants(w, maf = maf, ad = ad)
    # add name allele, seq.length, md5
    full_table <-
        ddply(ww, ~locus, function(x) {
            all_al <- unique(x$sequence)
            al_names <- {
                alexp <- ceiling(log10(length(all_al)))
                maxal <- 10^alexp
                numal <- seq_len(maxal)
                z <- str_pad(numal, alexp + 1, pad = "0")
                z[seq_along(all_al)]
            }
            # create df with alleles
            vdigest <- Vectorize(digest)
            df_alleles <-
                data.frame(
                    sequence = all_al,
                    variant = al_names,
                    nt = vapply(all_al, nchar, numeric(1)), # sequence length
                    md5 = vapply(all_al, vdigest, character(1))
                ) # md5 hash
            # join variant table with allele table
            left_join(x,
                df_alleles,
                by = "sequence"
            )
        }) |>
        select(
            .data$sample, .data$locus, .data$variant, .data$reads,
            .data$nt, .data$md5, .data$sequence
        ) |>
        arrange(.data$sample, .data$locus, .data$variant)
    full_table <- as_tibble(full_table)
    attr(full_table, "maf") <- attr(ww, "maf")
    attr(full_table, "ad") <- attr(ww, "ad")
    return(full_table)
}
