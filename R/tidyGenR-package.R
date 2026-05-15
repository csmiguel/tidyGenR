#' @keywords internal
"_PACKAGE"

## usethis namespace: start#
#' @import dplyr ggplot2 patchwork
#' @importFrom Biostrings DNAString DNAStringSet matchPattern readDNAStringSet
#'  width writeXStringSet reverseComplement
#' @importFrom dada2 dada derepFastq filterAndTrim getDadaOpt learnErrors
#'  loessErrfun makeSequenceTable mergePairs plotQualityProfile
#'  removeBimeraDenovo
#' @importFrom DECIPHER AlignSeqs BrowseSeqs AlignPairs
#' @importFrom digest digest
#' @importFrom glue glue
#' @importFrom methods is
#' @importFrom plyr daply ddply dlply ldply mapvalues
#' @importFrom readr write_delim
#' @importFrom ShortRead countFastq
#' @importFrom stats as.dist cmdscale loess predict reorder setNames
#' @importFrom stringr str_extract str_extract_all str_pad str_remove
#'  str_remove_all str_split str_which str_split_i
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column tibble
#' @importFrom tidyr as_tibble drop_na pivot_longer pivot_wider replace_na
#'  separate_longer_delim separate_wider_delim unite
#' @importFrom tidyselect all_of any_of everything where
#' @importFrom utils capture.output combn read.table write.table
#' @importFrom writexl write_xlsx
## usethis namespace: end
NULL
