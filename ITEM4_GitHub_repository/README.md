# ITEM4 — Seawater injection and electrokinetic sediment remediation

Reproducibility materials for the study **“Seawater injection and electrokinetic sediment remediation”**. This repository organizes the supplied R analysis script and preserves its historical figure-code appendix.

> **Status: partial reproducibility package.** The original sequencing files, environmental metadata file, SILVA reference databases, external Tax4Fun results, and PAST SIMPER output were not included in the provided upload. Consequently, the complete workflow has **not** been executed or validated end-to-end. Sections 10–11 remain pending verification.

## Repository contents

| Path | Description |
|---|---|
| `R/ITEM4_master_analysis.R` | Executable main workflow (sections 00–12, subject to required input files and checks below) |
| `archive/ITEM4_original_figure_blocks.R` | Commented historical figure-code blocks, for audit only |
| `archive/ITEM4_complete_original_upload.R` | Unmodified source script supplied by the authors |
| `input/` | Expected input folder structure; data not included |
| `output/` | Generated tables, figures, sequence files, and reproducibility metadata |
| `docs/REPRODUCIBILITY_NOTES.md` | Known limitations and required checks |
| `docs/DATA_AVAILABILITY_TEMPLATE.md` | Template for manuscript repository statement |

## Experimental groups

The source script describes **E-WI** (electrokinetic treatment with seawater injection), **WI** (seawater injection), and **Control**, with a common **T0** reference. Sampling positions H and O are distinguished in the treatments; the environmental dataset has 11 records (including T0), while the microbiome analysis excludes sample `V1OT1` after quality control, leaving nine 16S samples.

## Analysis workflow

1. Paired-end FASTQ filtering, denoising, merging and chimera removal using DADA2.
2. Taxonomic assignment with SILVA 138.1 and neighbor-joining phylogeny.
3. Phyloseq object assembly and sample QC.
4. Physicochemical (pH, ORP, C/H/N/S) and enzymatic analyses (MCA, MUF).
5. Alpha diversity (Observed ASVs, Shannon, Simpson).
6. Weighted and unweighted UniFrac PCoA.
7. Order-level taxonomic composition.
8. Order-level abundance heatmap.
9. Bray–Curtis nMDS, envfit and univariate dbRDA.
10. Mantel/Spearman taxon–environment associations (**pending verified implementation**).
11. Tax4Fun functional prediction and PAST SIMPER (**external results; not rerun here**).
12. Session information and package citation export.

## Input files required

Place the following files under `input/` before running:

```text
input/
├── dataset_env_2026.txt
├── fastq/
│   ├── <sample>_R1_trimmed.fastq[.gz]
│   └── <sample>_R2_trimmed.fastq[.gz]
└── reference/
    ├── silva_nr99_v138.1_train_set.fa.gz
    └── silva_species_assignment_v138.1.fa.gz
```

The script expects paired FASTQ filenames to share the same sample identifier before the first underscore. Metadata sample IDs must match FASTQ sample IDs. The metadata must include the variables expected by the analysis, including `treatment`, `electrode`, `Samples`, `pH`, `ORP`, `N_mean`, `C_mean`, `H_mean`, `S_mean`, `MCA_mean`, and `MUF_mean` where applicable. Check the exact column names and contents before execution.

**Do not upload raw sequencing reads, identifying data, or third-party reference databases unless their distribution terms permit it.** Link raw reads to their official public archive accession when available.

## Running the analysis

Install R and the packages below (Bioconductor packages require BiocManager):

```r
install.packages(c('BiocManager', 'ggplot2', 'dplyr', 'tidyr', 'tibble',
                   'patchwork', 'ggrepel', 'ape', 'vegan', 'phangorn'))
BiocManager::install(c('dada2', 'phyloseq', 'msa'))
```

From the repository root:

```bash
Rscript R/ITEM4_master_analysis.R
```

The workflow writes output to `output/figures/`, `output/tables/` and `output/`. Figures are exported as PDF and SVG. The script does not automatically install packages. Run it only after checking inputs and package compatibility.

## Reproducibility and interpretation

- The historical neighbor-joining tree is retained rather than substituting a model-optimized tree.
- T0 is a **single shared reference**, not independent biological replication across treatments.
- Reported SD values should not be interpreted as treatment-level biological replication.
- Negative blank-corrected enzymatic measurements are preserved.
- The source script describes SIMPER as performed in **PAST**, not R.
- The original archive contains commented historical alternatives; they are not all executable.

See [reproducibility notes](docs/REPRODUCIBILITY_NOTES.md) before citing the workflow as fully reproducible.

## Citation and data availability

Add the final article citation, public sequencing accession, and repository release/DOI once confirmed. See [data availability template](docs/DATA_AVAILABILITY_TEMPLATE.md).

## License

Code is released under the MIT License (see `LICENSE`). This license does not automatically apply to third-party reference datasets or unpublished experimental data.
