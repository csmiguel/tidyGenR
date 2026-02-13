#' Remove empty FASTQ
#'
#' Removes FASTA/FASTQ files with a number of reads less than or equal to
#' \code{min_reads} from a specified directory.
#'
#' @param path2fastq Character. Path to the directory containing FASTQ files.
#' @param pattern Character. Pattern used to identify FASTQ files in
#'   \code{path2fastq}.
#' @param min_reads Integer. FASTQ files with a number of reads less than or
#'   equal to this threshold are removed.
#'
#' @return
#' No return value. Called for its side effects: FASTA/FASTQ files in
#' \code{path2fastq} with a number of reads less than or equal to
#' \code{min_reads} are removed from disk. Invisibly returns \code{NULL}.
#' @examples
#' tmp <- tempdir()
#' fq <- file.path(tmp, "no_reads.fastq")
#' file.create(fq)
#' remove_poor_fastq(tmp, pattern = "fastq", min_reads = 0)
#' @export
remove_poor_fastq <- function(
    path2fastq, pattern = "fastq.gz",
    min_reads = 0) {
    fastq <- sort(list.files(
        path = path2fastq,
        pattern = pattern,
        full.names = TRUE
    ))
    no_reads <-
        vapply(fastq, function(x) {
            countFastq(x)$record
        }, numeric(1)) <= min_reads
    file.remove(fastq[no_reads])
    message(
        sum(!no_reads),
        " files have been KEPT.\n",
        sum(no_reads),
        " files REMOVED:\n",
        paste(fastq[no_reads], collapse = "\n")
    )
    invisible(NULL)
}
