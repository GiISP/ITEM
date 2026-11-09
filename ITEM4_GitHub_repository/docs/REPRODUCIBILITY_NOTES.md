# Reproducibility audit and outstanding tasks

## Verified from the supplied script (not from execution)

- Main script sections 00–12 were extracted from the 6,942-line upload without rewriting scientific calculations.
- Historical appendix is preserved separately and the complete upload is archived verbatim.
- The environmental values are embedded in the main script, including the common T0.
- The original code expects 10 samples before exclusion and nine afterward; these are expectations, not independently verified outputs.
- The script specifies SILVA 138.1, DADA2 pseudo-pooling, a neighbor-joining tree, and UniFrac/Bray–Curtis analyses.

## Before claiming full reproducibility

1. Supply and validate the FASTQ files, `dataset_env_2026.txt`, and both SILVA files.
2. Verify metadata field names and types, especially `treatment`, `electrode`, `Samples` and the eight environmental predictors.
3. Check sample identifiers and read-tracking QC; the source notes historical ASV dimensions only as targets.
4. Review the `envfit` result export: `env_fit$vectors$r` is labelled `r2` in the supplied code; check whether squaring is required for the intended reported statistic.
5. Verify environmental data consistency against the independent sequencing metadata audit output.
6. Restore and validate the original Mantel/Spearman test blocks and the FDR correction families; the supplied executable section currently only emits a pending message.
7. Provide external Tax4Fun functional prediction tables and original PAST SIMPER output; no equivalent R analysis has been substituted.
8. Review journal-specific figures, color palettes, taxon collapsing and axis labels before publication.
9. Run the full pipeline and inspect generated figures, logs, package versions and `sessionInfo.txt`.
10. Add the correct public sequence accession, authorship, manuscript reference, repository URL and DOI.

## Known constraints

- No raw FASTQ, metadata, reference FASTA or external prediction tables accompanied the uploaded script.
- No R runtime is available in the current packaging environment, so the R script has not been executed here.
- The original code uses relative paths: run it from the repository root.
- Figure outputs and tables are generated only when the required steps run successfully.
- The source appendix is preserved as commented archival code and should not be treated as a tested script.
