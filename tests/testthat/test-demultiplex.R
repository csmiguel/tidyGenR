#data("primers")

#### ATENTION ####
# this test can only be run locally. I have commented all code to prevent it
# from running in Github actions. The problem arises from when the fact that
# the the test will create snapshots where local paths are written. These
# paths cause a conflict when evaluating the snapshots in github actions.

# # download raw data
# raw <- system.file("extdata", "raw", package = "tidyGenR")
#
# # demultipex by locus single end data
# test_that("cutadapt script sh for single-end data is well formed", {
#     sh <- file.path(tempdir(), "cutadapt_se.sh")
#     suppressWarnings(
#         demultiplex(
#             cutadapt = "cutadapt",
#             sh_out = sh,
#             freads = list.files(raw,
#                 pattern = "1.fastq.gz",
#                 full.names = TRUE
#             ),
#             primers = primers[1:3, ],
#             mode = "se",
#             temp_folder = tempdir(),
#             outdir = file.path(tempdir(), "demultiplexed"),
#             run = FALSE
#         )
#     )
#     rl <- readLines(sh)
#     rlsub <- gsub(tempdir(), "any_path", rl)
#     out_rl <- file.path(tempdir(), "out_rl")
#     writeLines(rlsub, con = out_rl)
#     expect_snapshot_file(out_rl, "expected_cutadapt_se_sh.sh")
# })
#
# # paired-end mode forms valid script and demultiplexing with cutadapt produces the expected output
# test_that("cutadapt script sh for paired-end data is well formed", {
#     sh <- file.path(tempdir(), "cutadapt_pe_run.sh")
#     suppressWarnings(
#         demultiplex(
#             cutadapt = "cutadapt",
#             sh_out = sh,
#             freads = list.files(raw,
#                 pattern = "1.fastq.gz",
#                 full.names = TRUE
#             ),
#             rreads = list.files(raw,
#                 pattern = "2.fastq.gz",
#                 full.names = TRUE
#             ),
#             primers = primers[1:3, ],
#             mode = "pe",
#             temp_folder = tempdir(),
#             log_out = file.path(tempdir(), "cutadapt_pe_run.log"),
#             outdir = file.path(tempdir(), "demultiplexed_pe_run"),
#             run = FALSE
#         )
#     )
#     # sh produced is well formatted
#     rl <- readLines(sh)
#     rlsub <- gsub(tempdir(), "any_path", rl)
#     out_rl <- file.path(tempdir(), "out_rl")
#     writeLines(rlsub, con = out_rl)
#     expect_snapshot_file(out_rl, "expected_cutadapt_pe_run.sh")
# })
