#' Rename alleles in a dataframe based on reference alleles
#' @param df_alleles Dataframe with at least the following columns:
#' 'allele' 'sequence'.
#' @param ref_alleles Dataframe with at least the same columns as 'df_alleles',
#' or path multifasta FASTA with allele names in fasta headers.
#' @param replacements T/F, if TRUE, a dataframe with the replacements
#' correspondence is returned.
#' Optionally, a path to a multifasta with fasta headers being allele names.
#' @return Dataframe as 'df_alleles' with renamed alleles as in 'ref_alleles'.
#' @examples
#' data("genotypes")
#' ref_al <-
#'     system.file("extdata/reference_alleles.fasta", package = "tidyGenR")
#' rename_alleles(df_alleles = genotypes, ref_alleles = ref_al)
#' @export
rename_alleles <- function(
    df_alleles, ref_alleles,
    replacements = FALSE) {
    # if reference alleles are in fasta, convert to DF
    if ("character" %in% class(ref_alleles)) {
        bs <- readDNAStringSet(ref_alleles)
        ref_alleles <- data.frame(allele = names(bs), sequence = as.character(bs))
    }
    mandatory_col <- c("allele", "sequence")
    # check mandatory colnames
    if (
        !all(mandatory_col %in% names(df_alleles) &
            mandatory_col %in% names(ref_alleles))) {
        stop(paste0(
            "Input dataframes need to have the following",
            " variables: '",
            paste(mandatory_col, collapse = "', '"), "'"
        ))
    }
    # reduced df, one row per allele
    red_df_alleles <-
        distinct(select(df_alleles, .data$allele, .data$sequence))
    red_ref_alleles <-
        distinct(select(ref_alleles, .data$allele, .data$sequence))
    # rename alleles
    renamed_alleles <-
        left_join(red_df_alleles,
            red_ref_alleles,
            by = "sequence", suffix = c("", "_ref")
        ) |>
        mutate(allele = if_else(is.na(.data$allele_ref),
            .data$allele,
            .data$allele_ref
        )) |>
        select(-.data$allele_ref)
    # create new column in original df with renamed alleles as in reference
    df_alleles_renamed_old <-
        df_alleles |>
        left_join(renamed_alleles,
            by = "sequence",
            suffix = c("old", "")
        )
    # check allele names are unique within each locus. Warn if there is a risk of
    # renaming an allele as an already named allele for a given locus.
    check_names <-
        daply(df_alleles_renamed_old, ~locus, function(x) {
            lal <- length(unique(x$allele)) # number of different alleles
            lse <- length(unique(x$sequence)) # number of distinct sequences
            lal == lse
        })
    if (!all(check_names)) {
        loci_string <-
        paste(levels(as.factor(df_alleles_renamed_old$locus))[!check_names],
             collapse = ", ")
        stop("The number of allele names and sequences differ for loci: ",
            loci_string)
    }
    # sets col order as in original df
    df_alleles_renamed <-
        select(df_alleles_renamed_old, names(df_alleles))
    # df with replacements
    replacements_df <-
        select(
            df_alleles_renamed_old, .data$locus,
            .data$alleleold, .data$allele
        ) |>
        distinct() |>
        setNames(c("locus", "old_name", "new_name"))
    # return either replacements or
    if (replacements) {
        return(replacements_df)
    }
    if (!replacements) {
        return(df_alleles_renamed)
    }
}
