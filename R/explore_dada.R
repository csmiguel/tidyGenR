#' Explore variants output by DADA2 in the parameter space
#'
#' DADA2 is run with a set of desired parameters for a set of FASTQ files. The element
#' 'clustering' output by [dada2::dada()] containing all relevant statistics
#' from cluster formation are output in tidy format as the first element. These
#' tidy results are plotted in 4 different ways.
#' @details 'OMEGA_A', 'BAND_SIZE' and 'pool = T/F' parameters can have a strong
#' effect in the variants called by [dada2::dada()]. This function explores the
#' relation of read count, relative frequency of the variants and the
#' 'birth_pval' assigned to the clusters for any given values of starting
#' parameters. Critical parameters 'pool' and 'BAND_SIZE' can be specified as
#' arguments. Additional parameters can be set with [dada2::setDadaOpt()].
#' The 'log(-log(birth_pval))' proposed in Rosen et al. (2012) is computed and
#' represented in plots. For represenation purposes, a virtually infinite value
#' of 'log(-log(birth_pval)) = 10', is assigned by default to the first cluster
#' and to 'birth_pval = 0'. Variants in each, sample/locus combination are
#' ranked by abundance and plotted in the legend. Ranks >= 3 are
#' named as "3". For biallelic markers, a rank >= 3 implies a likely false
#' positive. This, visual representation can be used to decide to tune
#' [variant_call()].
#'  - _omega_a_: threshold for variants to be significant overabundant
#' 'log(-log(birth_pval))' (see Rosen et al. 2012). For exploration, it is
#' recommended to run [explore_dada()] with a large `omega_a`.
#'  - _band_size_: positive numbers set a band size in Needleman-Wunsch alignments.
#' In this context, ends free alignment is performed.
#' Zero turns off banding, triggering full Needleman-Wunsch alignments,
#' in which gapless alignment is performed
#' (see [issue](https://github.com/benjjneb/dada2/issues/1982)).
#'  - _pool_: calling variants pooling samples can increase sensitivity
#'  (see [dicussion](https://benjjneb.github.io/dada2/pseudo.html)).
#' @references
#' Rosen et al. (2012). _Denoising PCR-amplified metagenome data_.
#'  BMC Bioinformatics, 13(1).
#' @param fs Character vector with paths to FASTQ files.
#' @param sample_locus Regex expression with groups to extract sample (group 1)
#'  and loci (group 2) from "(^\[a-zA-Z0-9\]*)_(\[a-zA-Z0-9\]*)".
#' @param value_na Numeric to replace 'NA' or infinite values assigned to 'pval'
#'  or 'birth_pval' in clustering element from 'dada-class'.
#' @param reduced If TRUE, a reduced number of columns is returned.
#' If FALSE, all columns from from 'dada-class' clustering are returned.
#' @param omega_a "OMEGA_A" passed to [dada2::dada()].
#' @param band_size "BAND_SIZE" passed to [dada2::dada()].
#' @param pool Passed to [dada2::dada()].
#' @param vline Numeric x-intersection to annotate in plots p1 and p3.
#' @param hline_fr Numeric y-intersection to annotate 'frequency' in plots.
#' @param p_titles character(4) with plot names.
#' Passed to 'ggtitle()' in plots p1:p4.
#' @returns
#' List with tidy 'dada-class' clustering element and plots.
#'    1. tidy_dada: tidy 'clustering' element from 'dada-class'
#'        merged across loci.
#'    2. p1: plot 1, frequency of variants (sample x locus) against
#' 'log(-log(birth_pval))'.
#'    3. p2: plot 2, read count of variants against their frequency.
#'    4. p3: plot 3, p1 facetted by locus.
#'    5. p4: plot 4, p2 facetted by locus.
#' @examples
#' fq <-
#'  list.files(system.file("extdata", "truncated",
#'                         package = "tidyGenR"),
#'                         pattern = "F_filt.fastq.gz",
#'             full.names = TRUE)
#' explore_dada(fq,
#'     value_na = 10,
#'     reduced = TRUE,
#'     pool = FALSE,
#'     vline = 2,
#'     hline_fr = 0.1,
#'     omega_a = 0.9,
#'     band_size = 16
#' )
#' @export
explore_dada <- function(fs,
                         sample_locus = "(^[a-zA-Z0-9]*)_([a-zA-Z0-9]*)",
                         value_na = 10, reduced = TRUE,
                         omega_a = 0.9,
                         band_size = 16,
                         pool = FALSE,
                         vline = NULL, hline_fr = NULL,
                         p_titles = NULL) {
    # get loci from filenames
    loci <-
        unique(sort(str_extract(basename(fs),
            pattern = sample_locus,
            group = 2
        )))
    dd <- # tidy dada
        lapply(loci, function(lc) {
            # assert loci are unique in matched files:
            # match any alphanumerics upstream and donwstream from locus name in file
            #   string.
            fs2 <-
                fs[str_which(
                    basename(fs),
                    paste0("[A-Za-z0-9]*", lc, "[A-Za-z0-9]*")
                )]
            # assert only one locus is mathced
            if (length(unique(str_extract(basename(fs2),
                sample_locus,
                group = 2
            ))) > 1) {
                stop("More than one locus are matched.")
            }
            dd1 <-
                dada(fs2,
                    OMEGA_A = omega_a,
                    selfConsist = TRUE,
                    pool = pool,
                    BAND_SIZE = band_size
                ) |>
                dada2list(basename(fs2))
            # combine data from dataframes from elements in dada list.
            df <-
                lapply(seq_along(dd1), function(x) {
                    y <-
                        dd1[[x]]$clustering |>
                        mutate(
                            loglogp = log(-log(.data$birth_pval)),
                            # variant frequency
                            vfreq = .data$abundance / sum(.data$abundance),
                            # replace infinites by value_na
                            loglogp = if_else(is.infinite(.data$loglogp),
                                value_na, .data$loglogp
                            ),
                            sample = names(dd1)[x] |>
                                str_extract(sample_locus, group = 1),
                            locus = names(dd1)[x] |>
                                str_extract(sample_locus, group = 2)
                        ) |>
                        replace_na(list(loglogp = value_na)) |>
                        arrange(desc(.data$loglogp))
                    mutate(y, rank = seq_len(nrow(y)))
                }) |>
                # combine dataframes from different samples
                do.call(what = "rbind")
            if (reduced) {
                df <-
                    select(
                        df, .data$sample, .data$locus, .data$sequence,
                        .data$rank, .data$loglogp, .data$abundance, .data$vfreq
                    )
            }
            return(df)
        }) |>
        # combine dataframes from different loci
        do.call(what = "rbind") |>
        # recode variant rank. If more than 3 then tree.
        ddply(~ sample + locus, function(x) {
            y <- arrange(x, desc(.data$abundance))
            mutate(y,
                rank = seq_len(nrow(y)),
                rank2 = if_else(.data$rank > 2, 3, .data$rank)
            )
        }) |>
        as_tibble()
    # freq vs loglogp
    pomega <-
        dd |>
        ggplot(aes(
            x = .data$loglogp,
            y = .data$vfreq,
            color = as.character(.data$rank2)
        )) +
        geom_point() +
        geom_vline(
            xintercept = vline, linewidth = 1,
            linetype = "dotted", color = "grey"
        ) +
        geom_hline(
            yintercept = hline_fr, linewidth = 1,
            linetype = "dotted", color = "grey"
        ) +
        xlab("log(-log(birth_pval))") +
        ylab("frequency") +
        ggtitle(p_titles[1]) +
        scale_color_manual(
            name = "rank abundance",
            values = c(
                "1" = "green",
                "2" = "darkgreen",
                "3" = "red"
            ),
            labels = c("1", "2", ">2")
        ) +
        theme_classic() +
        theme(axis.text.x = element_text(
            angle = 60,
            hjust = 1,
            vjust = 1
        ))
    # freq vs abundance
    pabF <-
        dd |>
        ggplot(aes(
            x = .data$abundance,
            y = .data$vfreq,
            color = as.character(.data$rank2)
        )) +
        geom_point() +
        scale_x_continuous(
            transform = "log2",
            breaks = c(1, 10, 100, 1000, 10000)
        ) +
        geom_hline(
            yintercept = hline_fr, linewidth = 1,
            linetype = "dotted", color = "grey"
        ) +
        ylab("frequency") +
        ggtitle(p_titles[2]) +
        scale_color_manual(
            name = "rank abundance",
            values = c(
                "1" = "green",
                "2" = "darkgreen",
                "3" = "red"
            ),
            labels = c("1", "2", ">2")
        ) +
        theme_classic() +
        theme(axis.text.x = element_text(
            angle = 60,
            hjust = 1,
            vjust = 1
        ))

    # 2 split by locus
    # fr loglogp
    pomega_loc <-
        pomega +
        facet_wrap(~ .data$locus) +
        ggtitle(p_titles[3])

    p_ab_fr_loc <-
        pabF +
        facet_wrap(~ .data$locus) +
        ggtitle(p_titles[4])
    # list to return
    lp <- list(
        tidy_dada = dd,
        p1 = pomega,
        p2 = pabF,
        p3 = pomega_loc,
        p4 = p_ab_fr_loc
    )
    attr(lp, "band_size") <- band_size
    attr(lp, "pool") <- pool
    attr(lp, "omega_a") <- omega_a
    return(lp)
}
