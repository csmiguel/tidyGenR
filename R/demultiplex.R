#' Demultiplex reads by locus
#'
#' Sequencing reads are demultiplexed by locus into separate FASTQ files
#' (locus.sample.\[1|2\].fastq.gz) based on the locus-specific primer
#' sequences using `cutadapt`.
#'
#' @details Additional trimming of RC primers in 3' end of reads or adapter
#'  trimming requires additional runs of cutadapt in trimming mode with '-a'
#' argument. Default "-O 12" restricts partial matching of primers to full
#'  sequences in most cases. By default primers are match at any place in the
#' DNA sequence but this behaviour can be changed by using anchored primers.
#' @param interpreter Path to interpreter.
#' @param cutadapt Path to cutadapt executable.
#' @param freads Character vector with file paths to forward reads.
#' @param rreads Character vector with file paths to reverse reads.
#' @param primers Dataframe with primers
#' @param sh_out File name for cutadapt command.
#' @param outdir Directory to write demultiplexed FASTQ files.
#' Created if it does not exist.
#' @param mode "pe", paired-end; "se", single-end; or "linked" (beta, non-tested),
#' for linked primers.
#' @param log_out Path to write cutadapt log file.
#' @param temp_folder Directory to save temp files.
#' @param overlap --overlap (see cutadapt documentation).
#' @param e -e (see cutadapt documentation).
#' @param extraArgs Other arguments for cutadapt (eg "--max-n 0
#' --max-expected-errors 3 --minimum-length=20").
#' @param run T/F, whether to run the script (T) or just write it (F).
#' @rdname demultiplex
#' @returns Demultiplexed 'sample.locus.\[12\].fastq.gz'.
#' @examples
#' data("primers")
#' freads <-
#'  list.files(system.file("extdata", "raw",
#'                         package = "tidyGenR"),
#'                         pattern = "1.fastq.gz",
#'             full.names = TRUE)
#' rreads <-
#'  list.files(system.file("extdata", "raw",
#'                         package = "tidyGenR"),
#'                         pattern = "2.fastq.gz",
#'             full.names = TRUE)
#' demultiplex(
#'     cutadapt = "cutadapt",
#'     freads = freads,
#'     rreads = rreads,
#'     primers = primers,
#'     run = FALSE)
#'
#' @export
demultiplex <- function(interpreter = "/bin/bash",
                        cutadapt = "cutadapt",
                        freads,
                        rreads = NULL,
                        primers = NULL,
                        sh_out = "demultiplex.sh",
                        outdir = "demultiplexed",
                        log_out = "cutadapt.log",
                        temp_folder = tempdir(),
                        mode = "pe",
                        overlap = 15,
                        e = 0.15,
                        extraArgs = "",
                        run = TRUE) {
    # CHECKS ####
    # check raw
    message("Running 'check_raw_reads()':\n")
    meta <- check_raw_reads(freads = freads, rreads = rreads)
    if (!all(meta$checks)) stop("Some checks failed.")
    message("All checks on RAW data passed.\n\n")
    # check interpreter
    if (!file.exists(interpreter)) {
        warning("Path to interpreter does not exist.")
        # check path to demultiplexed exists or create
    }
    # sample names
    samples <- meta$samples
    # check primers
    message("Checking valid 'primers' dataframe.\n")
    check_primers(primers)
    message("All checks run.\n\n")
    # END checks #
    # START write samples + primers temp files ####
    if (!dir.exists(temp_folder)) {
        dir.create(temp_folder, recursive = TRUE)
    }
    # samples
    samples_file <- file.path(temp_folder, "samples")
    writeLines(samples, samples_file)
    # primers
    fw <- DNAStringSet(primers$fw)
    names(fw) <- primers$locus
    fwf <- file.path(temp_folder, "fw_primers")
    rv <- DNAStringSet(primers$rv)
    names(rv) <- primers$locus
    rvf <- file.path(temp_folder, "rv_primers")
    # END write files #

    if (mode != "linked") {
        writeXStringSet(fw, filepath = fwf)
        writeXStringSet(rv, filepath = rvf)
    } else if (mode == "linked") { # beta
        # write linked primers
        linked_primers <- as.character()
        for (i in seq_len(nrow(primers))) {
            linked_primers[i] <-
                paste0(
                    ">", primers$locus[i], "\n",
                    primers$fw[i], "...", primers$rv[i], "\n"
                )
        }
        writeLines(linked_primers, con = fwf, sep = "")
    }

    # create variables
    g <- paste0("file:", fwf)
    G <- paste0("file:", rvf)
    o <- file.path(outdir, "$sample.{name}.1.fastq.gz")
    p <- file.path(outdir, "$sample.{name}.2.fastq.gz")
    extf <- gsub("^[a-zA-Z0-9]*", "", meta[["f_reads"]][1])
    extr <- gsub("^[a-zA-Z0-9]*", "", meta[["r_reads"]][1])
    readsfw <- file.path(dirname(freads)[1], paste0("\"$sample\"", extf))
    if (!is.null(rreads)) {
        readsrv <- file.path(dirname(rreads)[1], paste0("\"$sample\"", extr))
    }
    # create cutdapt command
    cmd <- cutadapt_command(
        mode = mode, interpreter = interpreter,
        samples_file = samples_file,
        cutadapt = cutadapt, extraArgs = extraArgs,
        overlap = overlap,
        e = e, g = g, G = G, o = o, p = p,
        readsfw = readsfw, readsrv = readsrv,
        log_out = log_out
    )
    # Write the command to the file
    writeLines(cmd, con = sh_out)
    message("Cutadapt script written to:")
    message(sh_out, "\n")
    # run script
    if (run) {
        # check cutadapt version
        cutadapt_v <- check_cutadapt_version(cutadapt)
        if (is.null(cutadapt_v))
            stop("'cutadapt' has not been found. Try setting 'run = F'")
        # create ourdir
        stopifnot(is.character(outdir) && nchar(outdir) > 0)
        if (!dir.exists(outdir)) {
            dir.create(outdir, recursive = TRUE)
        }
        # make executable
        system2("chmod", args = c("+x", sh_out))
        message("Running cutadapt...")
        system(normalizePath(sh_out))
        message("Demultiplexing using primers from:")
        message(c(
            length(primers$locus),
            "loci: ",
            paste(primers$locus,
                collapse = " "
            )
        ))
        message("\nSamples demultiplexed:")
        message(c(
            length(samples),
            " samples: ",
            paste(samples, collapse = " ")
        ))
        message(
            "\nInto individual fastq.gz ",
            "files 'sample.locus.[1|2].fastq.gz' and written to:"
        )
        message(outdir, "\n")
    }
}

#' Create cutadapt command
#' @param samples_file Text file with samples names.
#' @param g Fasta file with forward primers.
#' @param G Fasta file with reverse primers.
#' @param o Path to demultiplexed forward reads.
#' @param p Path to demultiplexed reverse reads.
#' @param readsfw Path to input forward reads.
#' @param readsrv Path to input reverse reads.
#' @return cutadapt command based on the 'mode" ('pe', 'se', 'linked').
cutadapt_command <- function(mode,
                             interpreter,
                             samples_file, cutadapt,
                             extraArgs, overlap, e,
                             g, G, o, p, readsfw,
                             readsrv, log_out) {
    common_cmd <- c(
        paste0("#!", interpreter),
        paste("cat", samples_file, "| while read sample"),
        "do",
        paste(cutadapt, "\\")
    )

    # Use switch to handle different cases
    specific_cmd <- switch(mode,
        "linked" = {
            c(
                paste("--discard-untrimmed --no-indels", extraArgs, "\\"),
                paste("--overlap", overlap, "\\"),
                paste("-e", e, "\\"),
                paste("-g", g, "\\"),
                paste("-o", o, "\\"),
                readsfw
            )
        },
        "pe" = {
            c(
                paste(
                    "--discard-untrimmed --pair-adapters --no-indels",
                    extraArgs, "\\"
                ),
                paste("--overlap", overlap, "\\"),
                paste("-e", e, "\\"),
                paste("-g", g, "\\"),
                paste("-G", G, "\\"),
                paste("-o", o, "\\"),
                paste("-p", p, "\\"),
                paste(readsfw, "\\"),
                readsrv
            )
        },
        "se" = {
            c(
                paste("--discard-untrimmed --no-indels", extraArgs, "\\"),
                paste("--overlap", overlap, "\\"),
                paste("-e", e, "\\"),
                paste("-g", g, "\\"),
                paste("-o", o, "\\"),
                readsfw
            )
        },
        stop("Invalid sequence type")
    )

    final_cmd <- c(
        common_cmd,
        specific_cmd,
        paste("done", ifelse(is.character(log_out), paste0("> ", log_out), ""))
    )

    return(final_cmd)
}

# check cutadapt version
check_cutadapt_version <- function(p_cutadapt) {
    if (file.exists(p_cutadapt)) {
        cutadapt_v <-
            system2(p_cutadapt, "--version", stdout = TRUE) |>
            as.numeric()
        if (cutadapt_v < 2.0) {
            warning("Install cutadapt >2.0")
        }
    } else if (!file.exists(p_cutadapt)) {
        cutadapt_v <- NULL
        warning("**cutadapt** has not been found in '", p_cutadapt, "'.")
    }
    return(cutadapt_v)
}