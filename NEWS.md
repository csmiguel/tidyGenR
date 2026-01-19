# tidyGenR News
## tidyGenR 0.1.4 (2026-01-16)
- 'demultiplex(mode = 'linked')' has been tested with real data. The R primers needed to be reverse-complemented. Fixed.
## tidyGenR 0.1.3 (2025-08-26)
- updated 'gen_wide2genalex()'. 'tidyr::separate_wider_delim()' replaced with base R code.
- 'ploidy' attribute is kept for remove_hemizygotes

## tidyGenR 0.1.2 (2025-08-26)
- updated 'gen_wide2structure()'. 'tidyr::separate_longer_delim()' could not handle correctly NA's and the code was replaced with base R.
- 'ploidy' attribute is kept for genotype conversions.

## tidyGenR 0.1.1 (2025-08-08)
- Update NEWs and README.

## tidyGenR 0.1.0 (2025-04-06)

### Initial Release
- Core functions available.
- Clean .git. 
- Stable version.
- Old .git was removed and fresh .git was initiated after checking the rcmdcheck did not return any error.

