# tidyGenR News

## tidyGenR 0.1.8 (2026-05-15)
- fixed compare_calls(), did not allow creads = T in "genotypes".
- new function align_variants_ref() runs pairwise alignments with local reference sequences and attaches alignment stats to tidy variants. Useful to spot homologous/non-homologous variants.

## tidyGenR 0.1.7 (2026-02-13)
- fixed CRAN checks: tempdir and tempfiles for examples and tests, remove examples for internal `:::` functions, add `@value` where missing.
- remove unused data objects msa and qlist.

## tidyGenR 0.1.6 (2026-02-10)
- first submission to CRAN.
- trunc_amp(), includes 'filt_name' argument to add flexibility to output filtered names.
- x-axis scale has changed in 'explore_dada()'
- optimization vignette added.
- 'out_popart()' internal 'write_delim(quote)' changed to "none".
- replace absolute paths by temp paths in examples.

## tidyGenR 0.1.5 (2026-01-22)
- bug fixed in 'variant_call_dada()', dada2 is always converted to list to prevent errors associated to loci in one sample only.

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

