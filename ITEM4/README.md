# Scripts for DADA2 processing and sediment microbiome analysis

**Study:** *Seawater injection and electrokinetic sediment remediation*  
**Repository:** [GiISP/ITEM](https://github.com/GiISP/ITEM)

This repository contains R code and documentation for sediment physicochemical measurements and 16S rRNA gene microbiome analyses. Sections follow the order of the manuscript.

> **Reproducibility note:** The code excerpts below illustrate steps from the master script; they are not standalone commands. The full workflow has not been validated end-to-end because required FASTQ files, metadata and reference databases are not included. Some manuscript-specific analyses remain pending verification.

## Contents

- [01. DADA2 processing](#01-dada2-processing)
- [02. Taxonomic assignment and phylogeny](#02-taxonomic-assignment-and-phylogeny)
- [03. Phyloseq object and quality control](#03-phyloseq-object-and-quality-control)
- [04. Environmental and enzymatic measurements](#04-environmental-and-enzymatic-measurements)
- [05. Alpha diversity](#05-alpha-diversity)
- [06. Beta diversity — UniFrac PCoA](#06-beta-diversity--unifrac-pcoa)
- [07. Taxonomic composition](#07-taxonomic-composition)
- [08. Taxonomic heatmap](#08-taxonomic-heatmap)
- [09. Bray–Curtis nMDS, envfit and dbRDA](#09-braycurtis-nmds-envfit-and-dbrda)
- [10. Mantel and Spearman associations](#10-mantel-and-spearman-associations)
- [11. Tax4Fun and SIMPER](#11-tax4fun-and-simper)
- [12. Reproducibility and software versions](#12-reproducibility-and-software-versions)

## Load required packages

```r
required <- c("dada2", "phyloseq", "msa", "phangorn", "ape", "vegan",
              "ggplot2", "dplyr", "tidyr", "tibble", "patchwork", "ggrepel")
missing_pkgs <- setdiff(required, rownames(installed.packages()))
if (length(missing_pkgs)) stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
invisible(lapply(required, library, character.only = TRUE))
```

## Project structure

```text
ITEM/
├── README.md
├── LICENSE
├── R/
│   └── ITEM4_master_analysis.R
├── archive/
│   ├── ITEM4_complete_original_upload.R
│   └── ITEM4_original_figure_blocks.R
├── docs/
│   ├── REPRODUCIBILITY_NOTES.md
│   └── DATA_AVAILABILITY_TEMPLATE.md
├── input/
│   ├── fastq/
│   └── reference/
└── output/
    ├── figures/
    └── tables/
```

**Full script:** [`R/ITEM4_master_analysis.R`](R/ITEM4_master_analysis.R)  
**Historical figure scripts:** [`archive/ITEM4_original_figure_blocks.R`](archive/ITEM4_original_figure_blocks.R)

---

## 01. DADA2 processing

Quality filtering, error learning, ASV inference, paired-end merging and chimera removal.

```r
# Read paired-end FASTQ files
fnFs <- sort(list.files(fastq_dir, pattern = "_R1_trimmed.fastq(\\.gz)?$", full.names = TRUE))
fnRs <- sort(list.files(fastq_dir, pattern = "_R2_trimmed.fastq(\\.gz)?$", full.names = TRUE))

# Filter, learn errors and infer ASVs
filter_track <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs,
  maxN = 0, truncQ = 0, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
errF <- dada2::learnErrors(filtFs, nbases = 5e8, multithread = TRUE)
errR <- dada2::learnErrors(filtRs, nbases = 1e9, multithread = TRUE)

# Merge and remove chimeras (after denoising and dereplication)
mergers <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtab.nochim <- dada2::removeBimeraDenovo(
  dada2::makeSequenceTable(mergers), method = "consensus")
```

[Back to contents](#contents)

---

## 02. Taxonomic assignment and phylogeny

SILVA 138.1 classification and the neighbor-joining phylogeny retained in the source workflow.

```r
taxa <- dada2::assignTaxonomy(seqtab.nochim, silva_train, multithread = TRUE)
taxa <- dada2::addSpecies(taxa, silva_species)

mult <- msa::msa(seqs, method = "ClustalW", type = "dna")
treeNJ <- ape::NJ(phangorn::dist.ml(phangorn::as.phyDat(mult, type = "DNA")))
```

[Back to contents](#contents)

---

## 03. Phyloseq object and quality control

Build the analysis object and exclude sample `V1OT1` as documented in the source.

```r
prok_data_all <- phyloseq::phyloseq(
  phyloseq::otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  phyloseq::phy_tree(treeNJ),
  phyloseq::tax_table(taxa),
  phyloseq::sample_data(meta)
)
prok_data <- phyloseq::prune_samples(
  phyloseq::sample_names(prok_data_all) != "V1OT1", prok_data_all)
```

[Back to contents](#contents)

---

## 04. Environmental and enzymatic measurements

pH, ORP, C, H, N, S, Leu-MCA and MUF-glu; the shared T0 reference is not independent treatment replication.

```r
# Data table is defined in the full script (11 environmental records).
write.csv(env_full, "output/tables/environment_full_mean_SD.csv", row.names = FALSE)

# Original composite figure code is preserved in the archive.
# See archive/ITEM4_original_figure_blocks.R
```

[Back to contents](#contents)

---

## 05. Alpha diversity

Observed ASVs, Shannon and Simpson indices.

```r
alpha <- phyloseq::estimate_richness(
  prok_data, measures = c("Observed", "Shannon", "Simpson"))
write.csv(alpha, "output/tables/alpha_diversity.csv")
```

[Back to contents](#contents)

---

## 06. Beta diversity — UniFrac PCoA

Weighted and unweighted UniFrac distances and principal coordinate analysis.

```r
uw <- phyloseq::distance(prok_data, method = "unifrac", weighted = FALSE)
ww <- phyloseq::distance(prok_data, method = "unifrac", weighted = TRUE)
pcoa_u <- phyloseq::ordinate(prok_data, method = "PCoA", distance = uw)
pcoa_w <- phyloseq::ordinate(prok_data, method = "PCoA", distance = ww)
```

[Back to contents](#contents)

---

## 07. Taxonomic composition

Order-level relative abundance profiles.

```r
ps_rel <- phyloseq::transform_sample_counts(prok_data, function(x) x / sum(x))
ps_order <- phyloseq::tax_glom(ps_rel, taxrank = "Order", NArm = FALSE)
order_long <- phyloseq::psmelt(ps_order)
```

[Back to contents](#contents)

---

## 08. Taxonomic heatmap

Order-level abundance heatmap using the same aggregated phyloseq object.

```r
order_wide <- as(phyloseq::otu_table(ps_order), "matrix")
if (phyloseq::taxa_are_rows(ps_order)) order_wide <- t(order_wide)
# Taxonomic labels and heatmap construction are in the full script.
```

[Back to contents](#contents)

---

## 09. Bray–Curtis nMDS, envfit and dbRDA

Community ordination and environmental associations; dbRDA is fitted one variable at a time in the main script.

```r
bray <- vegan::vegdist(comm, method = "bray")
nmds <- vegan::metaMDS(bray, k = 2, trymax = 100, trace = FALSE)
env_fit <- vegan::envfit(nmds, env[, env_vars, drop = FALSE], permutations = 9999)
```

[Back to contents](#contents)

---

## 10. Mantel and Spearman associations

Historical code is preserved for scientific verification before use.

**Status:** Requires verification against the final manuscript analysis and original external outputs.

[Back to contents](#contents)

---

## 11. Tax4Fun and SIMPER

Functional predictions and SIMPER were performed outside the master R workflow; original outputs and methods need verification.

**Status:** Requires verification against the final manuscript analysis and original external outputs.

[Back to contents](#contents)

---

## 12. Reproducibility and software versions

Save R session information and package citations for the analysis environment.

```r
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
```

[Back to contents](#contents)

---

## Input files and execution

The master script expects trimmed paired-end FASTQ files in `input/fastq/`, SILVA 138.1 taxonomy reference files in `input/reference/`, and `input/dataset_env_2026.txt`. Supply these locally before running.

```r
source("R/ITEM4_master_analysis.R")
```

**Do not interpret this as an end-to-end tested command:** execution depends on missing inputs and verification of the source workflow.

## Data provenance and limitations

- The environmental table contains 11 records including the common initial T0.
- The microbiome workflow retains nine samples after excluding `V1OT1`.
- Historical figure blocks are archived, not standalone runnable scripts.
- Mantel/Spearman and Tax4Fun/SIMPER sections are flagged for verification; SIMPER was performed in PAST, not substituted with a new R analysis.
- See [`docs/REPRODUCIBILITY_NOTES.md`](docs/REPRODUCIBILITY_NOTES.md) for details.

## License

MIT License — see [`LICENSE`](LICENSE).
