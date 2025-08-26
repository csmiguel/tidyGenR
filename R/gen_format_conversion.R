#' Conversion among genotype formats: tidy, compact, wide, STRUCTURE
#'
#' @param gen Genotypes: 'tidy' in 'gen_tidy2compact()', 'gen_tidy2wide()',
#' 'gen_tidy2genalex()', and gen_tidy2integers()';
#' 'compact' in 'gen_compact2wide()'; 'wide' in 'gen_wide2structure()'.
#' @param write_out File to write structure formatted genotypes
#' (gen_wide2structure),
#' or GENALEX-formatted genotypes (in 'gen_tidy2genalex').
#' In 'gen_tidy2genalex', genotypes are written to tab-delimited
#' "txt" or "xlsx", depending on the extension ("txt", "xlsx").
#' @param popdata Dataframe with 'sample' and 'population' columns.
#' Populations are added to STRUCTURE output.
#'  If FALSE (Default), popdata is not added to STRUCTURE output.
#' Mandatory in 'gen_tidy2genalex'.
#' @param ploidy Genotypes ploidy in 'gen' (only for 'gen_tidy2genalex()').
#' @param delim Allele delimiter in genotype calls. Default "/". E.x. "AA/BB".
#' @param data_name Name of dataset to print to GENALEX xlsx.
#' @name genotype_conversion
#' @description
#' Conversion of genotype data between formats.
#' - *gen_tidy2compact*, from *tidy* to *compact*.
#' - *gen_compact2wide*, from *compact* to *wide*.
#' - *gen_tidy2wide*, from *tidy* into *wide*.
#' - *gen_wide2structure*, from *wide* to STRUCTURE-formatted text file.
#' - *gen_tidy2genalex*, from *tidy* to GENALEX-formatted *xlsx* or text file.
#' - *gen_tidy2integers*, recode 'allele' in *tidy* genotypes as integers.
#' @details
#' Genotypes from genotype() are returned as *tidy* data. Tidy data implies
#' one data point per row. Each row from *tidy* genotypes represents an allele
#' call for a given sample and locus. Rows with missing data are excluded.
#' Thus, *tidy* genotypes contain at least the three columns: *sample*,
#' *locus* and *allele*; bit often contain more columns (eg, *reads*,
#' *allele_no*, *nt*, *sequence*, *md5*, *population*, etc.). The *tidy*
#' format can be expanded column-wise and row-wise without altering the data
#' structure, and it can be easily manipulated. Some of the more handy functions
#'  here-included involve genotype conversion between
#' different formats:
#'  - *compact*, each row correspond to the genotype call for a given
#'  locus in a sample.
#'  Therefore, it only contains three columns, *sample*, *locus* and
#'  *genotype*. Default
#'  separator for alleles under *genotype* is '/'. For diploid data, hemizygotes
#'  are recoded as "allele1" instead of "allele1/NA".
#'  All other column metadata is lost.
#'  - *wide*, genotypes are recoded in samples (rows) x loci (columns)
#'  dataframe.
#'  The first column *samples*, contain sample names. All other columns contain
#'  loci names.
#'  Each cell contains genotypes in the *A/B* format.
#'
#' For the *STRUCTURE* format, '-9' are introduced as missing data in
#' STRUCTURE output.
#' Ploidy is detected automatically from 'allele_no'.
#' NAs are introduced if no allele_no found.
#' POPULATION of origin of the sample or other grouping factor such as SPECIES,
#'  are treated as follows:
#' @return
#' - *gen_tidy2compact*, *compact* genotypes with at least 'sample' ,'locus'
#'  and 'genotype' columns.
#' - *gen_compact2wide*, *wide* genotypes.
#' - *gen_tidy2wide*, *wide* genotypes.
#' - *gen_wide2structure*, plain text file formatted for STRUCTURE.
#' - *gen_tidy2genalex*, 'xls' file with genotype data for GENALEX.
#' - *gen_tidy2integers*, tidy genotypes with alleles recoded as
#' integers (1..n).
#' @examples
#' data("genotypes")
#' gencom <- gen_tidy2compact(genotypes)
#' gen_compact2wide(gencom)
#' gen_tidy2wide(genotypes)
#' gen_tidy2integers(genotypes)
#' # read metadata for sample populations
#' meta <-
#'     read.csv(system.file("extdata/metadata.csv", package = "tidyGenR"))
#' gen_wide2structure(gen_tidy2wide(genotypes),
#'     write_out = "genotypes.str", popdata = meta
#' )
#' gen_tidy2genalex(genotypes,
#'     popdata = meta,
#'     write_out = "dataset1.txt",
#' )
#' @export
gen_tidy2compact <- function(gen, delim = "/") {
    mand_vars <- c("locus", "sample", "allele")
    if (!all(mand_vars %in% names(gen))) {
        stop(
            "column names of 'gen' must contain ",
            paste(mand_vars, collapse = " ")
        )
    }
    if (!"allele_no" %in% names(gen)) {
        gen <- add_allele_no(gen)
    }
    z <-
        pivot_wider(gen,
            id_cols = c("sample", "locus"),
            names_from = "allele_no",
            values_from = "allele"
        )
    w <- unite(z, "genotype", -c(.data$sample, .data$locus),
        sep = delim, na.rm = TRUE
    )
    attr(w, "ploidy") <- attr(gen, "ploidy")
    return(w)
}
#' @rdname genotype_conversion
#' @export
gen_compact2wide <- function(gen) {
    mand_vars <- c("locus", "sample", "genotype")
    if (!all(mand_vars %in% names(gen))) {
        stop(
            "column names of 'gen' must contain: ",
            paste(mand_vars, collapse = ", ")
        )
    }
    w <-
      pivot_wider(gen,
        id_cols = "sample",
        names_from = "locus",
        values_from = "genotype"
    )
    attr(w, "ploidy") <- attr(gen, "ploidy")
    return(w)
}

#' @rdname genotype_conversion
#' @export
gen_tidy2wide <- function(gen) {
    gencomp <- gen_tidy2compact(gen)
    gen_compact2wide(gencomp)
}
#' @rdname genotype_conversion
#' @export
gen_wide2structure <- function(gen, write_out, popdata = FALSE, delim = "/") {
  if(isTRUE(attr(gen, "ploidy") == 2)) {
    message("All loci are assumed to have a ploidy of 2")
  } else {
    stop("The genotypes must have a 'ploidy' of 2 in their attributes: attr(gen, 'ploidy') == 2.")
  }
    gen <- column_to_rownames(gen, "sample")
    # separate alleles into two contigous rows
    h <-
      plyr::alply(gen, 1, function(x) {
        data.frame(a = str_split_i(x, "//?", 1),
                   b = str_split_i(x, "//?", 2)) |>
          t() |>
          as.data.frame()
      }) |>
      bind_rows() |>
      mutate(across(everything(), ~as.integer(as.factor(.x)))) |>
      setNames(names(gen))
    h[is.na(h)] <- "-9"
    # duplicate sample names
    snames <- rep(rownames(gen), each = 2)
    # dataframe with sample + pop
    if (inherits(popdata, "data.frame")) {
      test_1 <- all(c("sample", "population") %in% names(popdata))
      if (!test_1) {
        stop(
          "'popdata' must have a column named 'sample' matching samples",
          " in 'gen', and a 'population' descriptor if POPDATA flag is used ",
          "in STRUCTURE")
      } else {
        message("popdata is formatted correctly.")
      }
      # reorder pops as per sample names
        df1 <-
            data.frame(
                sample = snames,
                pops = popdata[["population"]][match(snames, popdata[["sample"]])]
            )
        df1[["pops"]] <- as.integer(as.factor(df1[["pops"]]))
        if (sum(is.na(df1$pops) > 0)) {
            stop("One or more samples did not match any of the populations.")
        }
    } else if (isFALSE(popdata)) {
      message("No popdata available.")
      df1 <-
        data.frame(sample = snames)
    } else {
      stop("'popdata' must be a dataframe.")
    }
    df2 <- cbind(df1, h)
    # first line in structure file
    a <- paste(colnames(df2)[which(colnames(df2) == colnames(h)[1]):ncol(df2)],
        collapse = " "
    )
    # lines with genotypes in structure file
    b <- apply(df2, 1, paste, collapse = " ")
    # remove if exists
    if (file.exists(write_out)) {
        file.remove(write_out)
    }
    # write to file
    Con <- file(write_out, open = "w")
    writeLines(a, con = Con)
    writeLines(b, con = Con)
    close(Con)
    message("STRUCTURE file written to ", write_out)
}

#' @rdname genotype_conversion
#' @export
gen_tidy2integers <- function(gen) {
    mand_vars <- c("sample", "locus", "allele")
    if (!all(mand_vars %in% names(gen))) {
        stop(
            "The mandatory columns '",
            paste(mand_vars, collapse = " "),
            "' have not been detected in 'gen'. Ensure the input is ",
            "a tidy genotype dataframe."
        )
    }
    z <-
        ddply(gen, ~locus, function(x) {
            x$allele <- as.integer(as.factor(x$allele))
            return(x)
        })
    return(tibble(z))
}

#' @rdname genotype_conversion
#' @export
gen_tidy2genalex <- function(gen, ploidy = 2,
                             popdata = FALSE,
                             write_out,
                             data_name = "dataset1") {
    # Check if the detected ploidy matches the set ploidy
    stopifnot(
        guess_ploidy(gen) == ploidy,
        !isFALSE(popdata),
        "data.frame" %in% class(popdata),
        all(c("sample", "population") %in% names(popdata))
    )
    # remove hemizygotes (if not it will cause errors in GENALEX)
    # and recode allele names as intergers for GENALEX.
    v <- gen_tidy2wide(gen_tidy2integers(remove_hemizygotes(gen)))
    w <- column_to_rownames(v, "sample")
    # 2 cols per locus
    gen2col <-
        separate_wider_delim(w,
            cols = everything(),
            delim = "/", names_sep = "_",
            too_few = "align_start"
        )
    # replace NAs with '0'. Genalex codes missing data as '0' for codominant data
    # and as '-1' for binary data.
    gen2col[is.na(gen2col)] <- "0"
    # sample names
    snames <- rownames(w)
    # order pops according to samples
    df1 <- data.frame(
        sample = snames,
        pops = popdata[["population"]][match(snames, popdata[["sample"]])]
    )
    # bind sample, population data and genotypes
    df2 <- cbind(df1, gen2col) |> arrange(.data$pops)
    # A1, no. loci; B1, no. samples; C1:Cn, samples per pop
    first_row <-
        c(
            ncol(w), nrow(w), length(unique(df2$pops)),
            as.numeric(table(df2$pops))
        )
    first_row <- c(first_row, rep("", ncol(df2) - length(first_row)))
    # second row
    second_row <- c(data_name, "", "", unique(df2$pops))
    second_row <- c(second_row, rep("", length(first_row) - length(second_row)))
    # third row
    third_row <- c("sample", "population", c(rbind(names(w), "")))
    df20 <- rbind(first_row, second_row, third_row, df2)
    # save to file
    extension <- tools::file_ext(write_out)
    if (extension == "xlsx") {
        write_xlsx(df20,
            path = write_out,
            col_names = FALSE
        )
    } else if (extension == "txt") {
        write.table(df20,
            write_out,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE, sep = "\t"
        )
    }
    message("GENALEX genotypes have been written to ", write_out)
}
