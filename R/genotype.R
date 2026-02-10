#' Genotype samples from tidy variants
#'
#' @details For 'ploidy':
#' - '1', force haploid genotypes.
#' - '2', force diploid genotypes. If number of alleles is == 1, ADt criterium
#' is applied to differentiate between "homozygotes" and "hemizygotes".
#' - 'poly', all variants are assumed to be real alleles.
#' No ADt filter is applied.
#' @param variants Tidy variants.
#' @param ploidy '1', force haploid genotypes; '2', force diploid genotypes;
#' 'poly', all variants are assumed to be real alleles.
#' @param ADt Threshold of AD to discriminate hemi/homo-zygotes.
#' @examples
#' data("variants")
#' genotype(variants[1:100,],
#'     ADt = 10, ploidy = 2
#' )
#' @export
#' @return Tidy genotypes. Dataframe with tidy genotypes.
#' 'sample, 'locus', 'allele', 'allele_no', 'reads', 'nt', 'md5', 'sequence'.
#' *reads*, supporting allele copies are corrected: 'reads' from tidy variants
#' are divided between resulting alleles in homozygotes and hemizygotes.
genotype <- function(variants,
                     ploidy = 2,
                     ADt = 10) {
    # checks
    mand_vars <- c("locus", "sample", "reads", "variant")
    stopifnot(
        "data.frame" %in% class(variants),
        all(mand_vars %in% names(variants))
    )
    variants$variant <- as.character(variants$variant)
    # genotype
    g <- geno_direct(variants = variants, ploidy = ploidy, ADt = ADt)
    # sort and arrange columns
    mand_names <- c("sample", "locus", "allele", "allele_no")
    v <- match(mand_names, names(g))
    v2 <- which(!names(g) %in% mand_names)
    g <-
        g[, c(v, v2)] |>
        arrange(across(all_of(mand_names))) |>
        as_tibble()
    # message ADt
    if (ploidy != 2) {
        message("ADt is ignored when ploidy != 2.")
    }
    # add attributes
    attr(g, "ploidy") <- ploidy
    if (!is.null(attr(variants, "MAF"))) {
        attr(g, "MAF") <- attr(variants, "MAF")
    }
    if (!is.null(attr(variants, "AD"))) {
        attr(g, "AD") <- attr(variants, "AD")
    }
    return(g)
}

#' Direct genotyping
#'
#' Performs straight genotyping from variants based on ploidy, ADt.
#' @rdname genotype
geno_direct <- function(variants, ploidy, ADt) {
    geno <-
        ddply(variants, ~ locus + sample, .drop = TRUE, function(x) {
            # .drop = T will drop unobserved combinations
            rreads <- x[["reads"]]
            names(rreads) <- x[["variant"]]
            if (length(rreads) == 0) {
                z <- NA
            } else if (ploidy == "poly") {
                z <- paste0(names(rreads), collapse = "/")
            } else if (ploidy == 1 & length(rreads) > 1) {
                z <- NA
            } else if (ploidy == 1 & length(rreads) == 1) {
                z <- names(rreads)
            } else if (ploidy == 2) {
                if (length(rreads) > 2) {
                    z <- NA
                } else if (length(rreads) == 2) {
                    z <- paste(names(rreads)[1], names(rreads)[2], sep = "/")
                } else if (length(rreads) == 1 & rreads > ADt) {
                    z <- paste(names(rreads)[1], names(rreads)[1], sep = "/")
                } else if (length(rreads) == 1 & rreads <= ADt) {
                    z <- paste(names(rreads)[1], "NA", sep = "/")
                }
            }
            data.frame(
                "locus" = unique(as.character(x$locus)),
                "sample" = unique(as.character(x$sample)),
                "genotype" = z
            )
        })
    # reformat to tidy
    a <-
        separate_longer_delim(geno, genotype, "/") |>
        add_allele_no() |>
        rename(variant = genotype)

    b <-
        rename(
            left_join(a, variants,
                by = c("sample", "locus", "variant")
            ),
            allele = .data$variant
        ) |>
        drop_na() |> # divide reads between alleles
        ddply(~ locus + sample, function(x) {
            y <-
                as.data.frame(table(x$allele)) |>
                setNames(c("allele", "freq"))
            left_join(x, y, by = "allele") |>
                rowwise() |>
                mutate(reads = .data$reads / .data$freq) |>
                select(-.data$freq)
        })
    return(b)
}

#' add allele_no to genotypes
#' @param gen genotypes
add_allele_no <- function(gen) {
    stopifnot(!"allele_no" %in% names(gen))
    gen_alleleno <-
        ddply(gen, ~ sample + locus, function(x) {
            mutate(x, allele_no = seq_len(nrow(x)))
        }) |>
        as_tibble()
    attr(gen_alleleno, "ploidy") <- attr(gen, "ploidy")
    return(gen_alleleno)
}
