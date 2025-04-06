#' Remove empty FASTQ
#'
#' Remove FASTA files with no reads in a given directory.
#' It considers both FASTQ and FASTA.
#'
#' @param path2fastq Folder with FASTQ files.
#' @param pattern Pattern to match FASTQ files in folder.
#' @param min_reads FASTQ files with a <= number of reads are deleted.
#' @examples
#' # create empty FASTQ
#' fq <- "no_reads.fastq"
#' system2("touch", fq)
#' # remove empty fastq
#' remove_poor_fastq(".", pattern = "fastq")
#' @export
remove_poor_fastq <- function(
    path2fastq = NULL, pattern = "fastq.gz",
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
}
