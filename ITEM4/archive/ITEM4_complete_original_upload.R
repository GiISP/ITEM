# ITEM4 — MASTER ANALYSIS SCRIPT (MANUSCRIPT ORDER)
# Manuscript: Seawater injection and electrokinetic sediment remediation
# SIMPER was run in PAST, not in R.

# FIGURE / RESULTS NAVIGATION ------------------------------------------------
# IMPORTANT: FIGURE SLOTS A-G indicate order of manuscript microbiology results,
# ENV. Environmental chemistry and enzymes -> Section 04
# A. Alpha diversity (Observed ASVs, Shannon) -> Section 05
# B. PCoA weighted/unweighted UniFrac -> Section 06
# C. Taxonomic composition -> Section 07
# D. Taxonomic heatmap -> Section 08 
# E. Bray-Curtis nMDS / envfit / dbRDA -> Section 09
# F. Mantel / Spearman taxon-environment links -> Section 10 
# G. Tax4Fun KO/pathway prediction and PAST SIMPER -> Section 11 
# Export rules: all figure-producing sections should write PDF + SVG, and
# export the underlying processed table to output/tables.


# ENVIRONMENTAL DATA: all 11 samples; pH, ORP, C, H, N, S, MCA, MUF; means and SD.
# Run from the project root with input/fastq, input/reference and input/dataset_env_2026.txt.
# Original analysis excerpts are preserved at the end as non-executing reference.

# 00. PROJECT SETUP -----------------------------------------------------------
required <- c('dada2','phyloseq','msa','phangorn','ape','vegan',
              'ggplot2','dplyr','tidyr','tibble','patchwork','ggrepel')
missing_pkgs <- setdiff(required, rownames(installed.packages()))
if(length(missing_pkgs)) stop('Install required packages: ',paste(missing_pkgs,collapse=', '))
invisible(lapply(required, library, character.only=TRUE))
set.seed(2026) # for stochastic downstream ordinations/permutations only
project <- '.'
input_dir <- file.path(project,'input')
fastq_dir <- file.path(input_dir,'fastq')
ref_dir <- file.path(input_dir,'reference')
output_dir <- file.path(project,'output')
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(output_dir,'figures'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(output_dir,'tables'),recursive=TRUE,showWarnings=FALSE)

# 01. ORIGINAL DADA2 PROCESSING -----------------------------------------------
# Original inputs: *_R1_trimmed.fastq and *_R2_trimmed.fastq.
# The original script has a stray comma in filterAndTrim(); removed as syntax fix.
fnFs <- sort(list.files(fastq_dir,pattern='_R1_trimmed.fastq(\\.gz)?$',full.names=TRUE))
fnRs <- sort(list.files(fastq_dir,pattern='_R2_trimmed.fastq(\\.gz)?$',full.names=TRUE))
stopifnot(length(fnFs)>0, length(fnFs)==length(fnRs))
# Original extraction takes first underscore-delimited token. Validate uniqueness.
sample_names_raw <- vapply(strsplit(basename(fnFs),'_'),`[`,character(1),1)
stopifnot(!anyDuplicated(sample_names_raw))
stopifnot(identical(sample_names_raw,
  vapply(strsplit(basename(fnRs),'_'),`[`,character(1),1)))
dir.create(file.path(output_dir,'filtered'),showWarnings=FALSE)
filtFs <- file.path(output_dir,'filtered',paste0(sample_names_raw,'_F_filt.fastq.gz'))
filtRs <- file.path(output_dir,'filtered',paste0(sample_names_raw,'_R_filt.fastq.gz'))
filter_track <- dada2::filterAndTrim(fnFs,filtFs,fnRs,filtRs,
  maxN=0,truncQ=0,rm.phix=TRUE,compress=TRUE,multithread=TRUE)
errF <- dada2::learnErrors(filtFs,nbases=5e8,multithread=TRUE,randomize=TRUE)
errR <- dada2::learnErrors(filtRs,nbases=1e9,multithread=TRUE,randomize=TRUE)
derepFs <- dada2::derepFastq(filtFs,verbose=TRUE)
derepRs <- dada2::derepFastq(filtRs,verbose=TRUE)
names(derepFs) <- names(derepRs) <- sample_names_raw
dadaFs <- dada2::dada(derepFs,err=errF,pool='pseudo',multithread=TRUE)
dadaRs <- dada2::dada(derepRs,err=errR,pool='pseudo',multithread=TRUE)
mergers <- dada2::mergePairs(dadaFs,derepFs,dadaRs,derepRs,verbose=TRUE)
seqtab <- dada2::makeSequenceTable(mergers)
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab,method='consensus',
  multithread=TRUE,verbose=TRUE)
getN <- function(x) sum(dada2::getUniques(x))
track <- cbind(filter_track,sapply(dadaFs,getN),sapply(dadaRs,getN),
               sapply(mergers,getN),rowSums(seqtab.nochim))
colnames(track) <- c('input','filtered','denoisedF','denoisedR','merged','nonchim')
rownames(track) <- sample_names_raw
write.csv(track,file.path(output_dir,'tables','DADA2_read_tracking.csv'))
# Historical checkpoints reported in original console log: 10 x 5029 ASVs before
# chimera removal, 10 x 1497 after. Treat as QC targets, not hardcoded results.
message('Sequence table dimensions: ',paste(dim(seqtab.nochim),collapse=' x '))

# 02. ORIGINAL SILVA ASSIGNMENT AND PHYLOGENY ---------------------------------
silva_train <- file.path(ref_dir,'silva_nr99_v138.1_train_set.fa.gz')
silva_species <- file.path(ref_dir,'silva_species_assignment_v138.1.fa.gz')
stopifnot(file.exists(silva_train),file.exists(silva_species))
taxa <- dada2::assignTaxonomy(seqtab.nochim,silva_train,multithread=TRUE)
taxa <- dada2::addSpecies(taxa,silva_species)
write.csv(taxa,file.path(output_dir,'tables','ASV_taxonomy_SILVA138_1.csv'))
dada2::uniquesToFasta(seqtab.nochim,
  fout=file.path(output_dir,'ASV_representative_sequences.fasta'),
  ids=colnames(seqtab.nochim))
seqs <- dada2::getSequences(seqtab.nochim)
names(seqs) <- seqs
mult <- msa::msa(seqs,method='ClustalW',type='dna',order='input')
phang.align <- phangorn::as.phyDat(mult,type='DNA')
dm <- phangorn::dist.ml(phang.align)
treeNJ <- ape::NJ(dm)
# Important provenance note: the original script optimized a GTR model in
# 'fitGTR' but built phyloseq with treeNJ, not fitGTR$tree.
# Therefore the original NJ tree is intentionally retained here.
ape::write.tree(treeNJ,file.path(output_dir,'ASV_tree_original_NJ.nwk'))

# 03. PHYLOSEQ OBJECT / QC ----------------------------------------------------
metadata_file <- file.path(input_dir,'dataset_env_2026.txt')
stopifnot(file.exists(metadata_file))
meta <- read.delim(metadata_file,header=TRUE,row.names=1,check.names=FALSE)
stopifnot(all(rownames(seqtab.nochim) %in% rownames(meta)))
meta <- meta[rownames(seqtab.nochim),,drop=FALSE]
prok_data_all <- phyloseq::phyloseq(
  phyloseq::otu_table(seqtab.nochim,taxa_are_rows=FALSE),
  phyloseq::phy_tree(treeNJ),phyloseq::tax_table(taxa),
  phyloseq::sample_data(meta))
# V1OT1 excluded in original analysis owing to insufficient reads after chimera filtering.
stopifnot('V1OT1' %in% phyloseq::sample_names(prok_data_all))
prok_data <- phyloseq::prune_samples(
  phyloseq::sample_names(prok_data_all)!='V1OT1',prok_data_all)
prok_data <- phyloseq::prune_taxa(phyloseq::taxa_sums(prok_data)>0,prok_data)
stopifnot(phyloseq::nsamples(prok_data)==9)
saveRDS(prok_data,file.path(output_dir,'phyloseq_9samples_original_tree.rds'))

# 04. ENVIRONMENTAL AND ENZYMATIC RESULTS (FIRST MANUSCRIPT FIGURES) ------
# Transcribed from the author-provided ENV_FULL mean +/- SD matrix.
# SDs describe analytical/within-mesocosm variation, NOT independent
# biological treatment replication. T0 is a common reference condition.
# V1OT1 is retained here although excluded from 16S analyses.
env_full <- tibble::tibble(
  Sample = c("T0", "V1HT1", "V1HT2", "V1OT1", "V1OT2", "V2HT1", "V2HT2", "V2OT1", "V2OT2", "V3CT1", "V3CT2"),
  Treatment = c("Initial", "E-WI", "E-WI", "E-WI", "E-WI", "WI", "WI", "WI", "WI", "Control", "Control"),
  Compartment = c("Initial", "H", "H", "O", "O", "H", "H", "O", "O", "C", "C"),
  Time = c("T0", "T1", "T2", "T1", "T2", "T1", "T2", "T1", "T2", "T1", "T2"),
  pH = c(7.88, 7.185, 8.01, 7.225, 8.045, 7.755, 8.005, 7.715, 7.98, 7.89, 8.025),
  pH_SD = c(NA_real_, 0.077782, 0.014142, 0.007071, 0.021213, 0.007071, 0.021213, 0.035355, 0.007071, 0.069706, 0.014142),
  ORP = c(-304.0, 39.0, 172.5, 231.5, 195.0, 77.0, 203.5, 46.5, 220.0, -238.5, -302.0),
  ORP_SD = c(NA_real_, 4.242641, 3.535534, 24.748737, 8.485281, 15.556349, 6.363961, 26.162951, 15.556349, 25.455844, 26.162951),
  N_mean = c(1.028723, 3.020518, 0.210098, 2.579916, 0.562667, 3.015412, 2.601548, 2.518925, 3.199029, 2.882339, 1.973543),
  N_SD = c(0.198926, 0.955795, 0.051387, 0.621599, 0.205611, 0.113501, 0.819283, 0.104063, 0.434538, 0.764624, 0.099705),
  C_mean = c(28.474261, 23.050563, 25.080316, 25.414129, 26.686532, 28.222157, 28.196222, 16.518648, 30.193817, 13.630208, 10.185317),
  C_SD = c(7.907091, 14.176583, 3.490981, 18.310202, 10.119413, 5.354117, 18.140989, 13.375134, 17.112742, 10.505006, 5.115372),
  H_mean = c(6.870688, 10.601628, 6.472906, 8.579289, 6.267412, 7.190078, 8.919028, 5.445488, 7.691238, 9.940284, 5.581613),
  H_SD = c(0.679111, 0.557827, 1.512672, 5.185893, 0.678218, 0.945823, 3.571357, 0.680335, 0.445441, 1.821833, 0.181358),
  S_mean = c(11.15311, 6.88161, 3.543277, 3.559849, 1.78608, 1.147717, 1.650637, 1.118691, 0.472783, 5.254055, 0.90899),
  S_SD = c(4.335021, 4.155086, 0.36287, 2.47933, 0.760169, 0.44156, 0.533341, 0.182447, 0.324619, 3.179452, 0.524806),
  MCA_mean = c(2205.620823, -0.815651, 1292.240888, 0.249755, 2419.061834, 501.33964, 10530.798717, 1071.033972, 11424.480669, 57.06027, 1420.022184),
  MCA_SD = c(1299.122448, 0.070701, 103.17213, 0.225873, 986.548443, 35.727846, 1212.514891, 139.771989, 1088.158038, 9.272893, 577.773614),
  MUF_mean = c(1069.441822, 1.880171, 222.662249, 5.856358, 331.573879, 1615.644592, 3602.388221, 1877.27906, 3292.032753, 303.654439, 1320.276313),
  MUF_SD = c(222.491194, 2.003459, 213.991801, 1.426911, 71.700262, 353.252914, 461.46729, 357.362212, 997.562247, 29.584931, 786.499993)
)
env_full$Time <- factor(env_full$Time, levels=c('T0','T1','T2'))
stopifnot(nrow(env_full)==11, !anyDuplicated(env_full$Sample))
write.csv(env_full,file.path(output_dir,'tables','environment_full_mean_SD.csv'),row.names=FALSE)

# One consistent figure per family; panel layout may need adaptation to the
# journal's final composite figures. Negative blank-corrected enzyme values
# are preserved; they are NOT truncated or log-transformed.
env_specs <- data.frame(
  variable=c('pH','ORP','N_mean','C_mean','H_mean','S_mean','MCA_mean','MUF_mean'),
  sd=c('pH_SD','ORP_SD','N_SD','C_SD','H_SD','S_SD','MCA_SD','MUF_SD'),
  label=c('pH','ORP (mV)','Nitrogen (% DW)','Carbon (% DW)',
          'Hydrogen (% DW)','Sulfur (% DW)',
          'Leucine aminopeptidase (nmol g-1 h-1)',
          'beta-D-glucosidase (nmol g-1 h-1)'))
plot_environment <- function(variable, sd_col, label) {
  x <- env_full
  x$group <- factor(paste(x$Treatment,x$Compartment,sep=' / '),
                    levels=unique(paste(x$Treatment,x$Compartment,sep=' / ')))
  x$sd_plot <- x[[sd_col]]
  x$mean_plot <- x[[variable]]
  ggplot2::ggplot(x,ggplot2::aes(x=group,y=mean_plot,color=Time)) +
    ggplot2::geom_point(position=ggplot2::position_dodge(width=.55),size=2.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_plot-sd_plot,
       ymax=mean_plot+sd_plot),width=.15,
       position=ggplot2::position_dodge(width=.55),na.rm=TRUE) +
    ggplot2::labs(x=NULL,y=label,color='Time') +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,hjust=1))
}
for (i in seq_len(nrow(env_specs))) {
  z <- env_specs[i,]
  g <- plot_environment(z$variable,z$sd,z$label)
  ggplot2::ggsave(file.path(output_dir,'figures',paste0('Environment_',z$variable,'.pdf')),
                  g,width=9,height=5)
  ggplot2::ggsave(file.path(output_dir,'figures',paste0('Environment_',z$variable,'.svg')),
                  g,width=9,height=5)
}
# The independent 16S metadata must agree with env_full on shared variables.
# Do not overwrite original sequencing metadata silently.
metadata_env_crosscheck <- merge(
  tibble::rownames_to_column(as.data.frame(meta),'Sample'),
  as.data.frame(env_full),by='Sample',suffixes=c('_16S','_ENV'))
write.csv(metadata_env_crosscheck,
  file.path(output_dir,'tables','environment_vs_16S_metadata_audit.csv'),row.names=FALSE)

# 05. MANUSCRIPT: MICROBIAL ALPHA DIVERSITY (FIGURE SLOT A) ---------------------------------------------------------
alpha <- phyloseq::estimate_richness(prok_data,measures=c('Observed','Shannon','Simpson'))
alpha$Sample <- rownames(alpha)
alpha <- dplyr::left_join(alpha,
  tibble::rownames_to_column(as.data.frame(phyloseq::sample_data(prok_data)),
                            var='Sample'),by='Sample')
write.csv(alpha,file.path(output_dir,'tables','alpha_diversity.csv'),row.names=FALSE)
alpha_long <- tidyr::pivot_longer(alpha,cols=c('Observed','Shannon'),
                                  names_to='Metric',values_to='Value')
g_alpha <- ggplot2::ggplot(alpha_long,ggplot2::aes(x=Sample,y=Value)) +
  ggplot2::geom_col(fill='grey45') + ggplot2::facet_wrap(~Metric,scales='free_y',ncol=1) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
ggplot2::ggsave(file.path(output_dir,'figures','Alpha_diversity.pdf'),g_alpha,width=9,height=7)
ggplot2::ggsave(file.path(output_dir,'figures','Alpha_diversity.svg'),g_alpha,width=9,height=7)

# 06. MANUSCRIPT: BETA DIVERSITY / UNIFRAC (FIGURE SLOT B) ------------------------------------------------
# Uses original NJ phylogenetic tree, with original count table.
uw <- phyloseq::distance(prok_data,method='unifrac',weighted=FALSE)
ww <- phyloseq::distance(prok_data,method='unifrac',weighted=TRUE)
write.csv(as.matrix(uw),file.path(output_dir,'tables','unweighted_unifrac.csv'))
write.csv(as.matrix(ww),file.path(output_dir,'tables','weighted_unifrac.csv'))
pcoa_w <- phyloseq::ordinate(prok_data,method='PCoA',distance=ww)
pcoa_u <- phyloseq::ordinate(prok_data,method='PCoA',distance=uw)
plot_pcoa <- function(ord,title) {
  phyloseq::plot_ordination(prok_data,ord,type='samples',color='treatment',
    shape='electrode',label='Samples') + ggplot2::labs(title=title) +
    ggplot2::theme_classic()
}
p_pcoa <- plot_pcoa(pcoa_w,'Weighted UniFrac') |
          plot_pcoa(pcoa_u,'Unweighted UniFrac')
ggplot2::ggsave(file.path(output_dir,'figures','PCoA_UniFrac.pdf'),
                p_pcoa,width=12,height=6)
ggplot2::ggsave(file.path(output_dir,'figures','PCoA_UniFrac.svg'),
                p_pcoa,width=12,height=6)

# 07. MANUSCRIPT: TAXONOMIC COMPOSITION (FIGURE SLOT C) ---------------------------------------------------
# Order-level relative abundances; no arbitrary ASV pre-filter before aggregation.
ps_rel <- phyloseq::transform_sample_counts(prok_data,function(x) x/sum(x))
ps_order <- phyloseq::tax_glom(ps_rel,taxrank='Order',NArm=FALSE)
order_long <- phyloseq::psmelt(ps_order)
write.csv(order_long,file.path(output_dir,'tables','order_relative_abundance_long.csv'),
          row.names=FALSE)
g_order <- ggplot2::ggplot(order_long,
   ggplot2::aes(x=Sample,y=Abundance,fill=Order)) +
   ggplot2::geom_col() + ggplot2::theme_classic() +
   ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1)) +
   ggplot2::labs(y='Relative abundance',x=NULL)
ggplot2::ggsave(file.path(output_dir,'figures','Taxonomic_orders_all.pdf'),
                g_order,width=11,height=6)
ggplot2::ggsave(file.path(output_dir,'figures','Taxonomic_orders_all.svg'),
                g_order,width=11,height=6)
# The manuscript-specific taxonomic palette and taxon-collapse rule must be
# copied from the final figure-producing block after visual audit.

# 08. MANUSCRIPT: TAXONOMIC HEATMAP (FIGURE SLOT D) ---------------------------
# Source: final order-level abundance table produced in Section 06.

# 08b. HEATMAP FROM ORIGINAL ORDER-LEVEL RELATIVE ABUNDANCES -----------------
# Uses the same order-level object as the taxonomic barplot. All taxa retained.
order_wide <- as(phyloseq::otu_table(ps_order), 'matrix')
if (phyloseq::taxa_are_rows(ps_order)) order_wide <- t(order_wide)
order_tax <- as.data.frame(phyloseq::tax_table(ps_order))
order_labels <- as.character(order_tax$Order)
order_labels[is.na(order_labels) | order_labels==''] <- 'Unclassified'
colnames(order_wide) <- make.unique(order_labels)
write.csv(order_wide, file.path(output_dir,'tables','order_heatmap_matrix.csv'))
heat_long <- as.data.frame(as.table(order_wide))
colnames(heat_long) <- c('Sample','Order','Abundance')
g_heat <- ggplot2::ggplot(heat_long,ggplot2::aes(x=Sample,y=Order,fill=Abundance)) +
  ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c() +
  ggplot2::theme_minimal(base_size=9) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
ggplot2::ggsave(file.path(output_dir,'figures','Fig05_Order_heatmap.pdf'),g_heat,width=10,height=12)
ggplot2::ggsave(file.path(output_dir,'figures','Fig05_Order_heatmap.svg'),g_heat,width=10,height=12)

# 09. MANUSCRIPT: BRAY-CURTIS / NMDS / ENVFIT / DBRDA (FIGURE SLOT E) ------------------------------------
comm <- as(phyloseq::otu_table(ps_rel),'matrix')
if(phyloseq::taxa_are_rows(ps_rel)) comm <- t(comm)
comm <- comm[phyloseq::sample_names(prok_data),,drop=FALSE]
env <- as.data.frame(phyloseq::sample_data(prok_data))
env <- env[rownames(comm),,drop=FALSE]
env_vars <- c('pH','ORP','N_mean','C_mean','H_mean','S_mean','MCA_mean','MUF_mean')
stopifnot(all(env_vars %in% colnames(env)))
env[env_vars] <- lapply(env[env_vars],function(x) as.numeric(gsub(',','.',as.character(x))))
stopifnot(!anyNA(env[env_vars]))
bray <- vegan::vegdist(comm,method='bray')
nmds <- vegan::metaMDS(bray,k=2,trymax=100,trace=FALSE)
env_fit <- vegan::envfit(nmds,env[,env_vars,drop=FALSE],permutations=9999)
write.csv(data.frame(variable=rownames(env_fit$vectors$arrows),
  r2=env_fit$vectors$r,p=env_fit$vectors$pvals),
  file.path(output_dir,'tables','envfit_results.csv'),row.names=FALSE)
# nMDS ordination with fitted environmental vectors.
nmds_sites <- as.data.frame(vegan::scores(nmds,display='sites'))
nmds_sites$Sample <- rownames(nmds_sites)
vec <- as.data.frame(vegan::scores(env_fit,display='vectors'))
vec$Variable <- rownames(vec)
write.csv(nmds_sites,file.path(output_dir,'tables','nMDS_sites.csv'),row.names=FALSE)
write.csv(vec,file.path(output_dir,'tables','nMDS_envfit_vectors.csv'),row.names=FALSE)
g_nmds <- ggplot2::ggplot(nmds_sites,ggplot2::aes(x=NMDS1,y=NMDS2,label=Sample)) +
  ggplot2::geom_point(size=3) + ggplot2::geom_text(vjust=-.7,size=3) +
  ggplot2::geom_segment(data=vec,ggplot2::aes(x=0,y=0,xend=NMDS1,yend=NMDS2),
                        inherit.aes=FALSE,arrow=grid::arrow(length=grid::unit(.15,'cm'))) +
  ggplot2::geom_text(data=vec,ggplot2::aes(x=NMDS1,y=NMDS2,label=Variable),
                     inherit.aes=FALSE,hjust=0,size=3) +
  ggplot2::coord_equal() + ggplot2::theme_classic() +
  ggplot2::labs(subtitle=paste('Bray-Curtis; stress =',round(nmds$stress,3)))
ggplot2::ggsave(file.path(output_dir,'figures','Fig06_nMDS_envfit.pdf'),g_nmds,width=9,height=7)
ggplot2::ggsave(file.path(output_dir,'figures','Fig06_nMDS_envfit.svg'),g_nmds,width=9,height=7)

# Univariate dbRDA only. Multivariable fitting is not justified with n=9
# without a separate, explicitly defined model selection strategy.
dbrda_one <- function(v) {
  f <- stats::as.formula(paste('comm ~',v))
  mod <- vegan::capscale(f,data=env,distance='bray',add=TRUE)
  a <- vegan::anova.cca(mod,permutations=9999)
  data.frame(variable=v,R2=unname(vegan::RsquareAdj(mod)$r.squared),
             p=unname(a$`Pr(>F)`[1]))
}
dbrda_table <- dplyr::bind_rows(lapply(env_vars,dbrda_one))
dbrda_table$FDR <- p.adjust(dbrda_table$p,method='BH')
write.csv(dbrda_table,file.path(output_dir,'tables','univariate_dbRDA.csv'),row.names=FALSE)

# 10. MANUSCRIPT: TAXON-ENVIRONMENT ASSOCIATIONS (FIGURE SLOT F) ------------------------------------------
# The original manuscript's Mantel/Spearman implementation and FDR families
# require an exact code-level comparison with the final analysis block.
# DO NOT silently replace those tests with simple cor() or newly chosen models.
message('ACTION REQUIRED: port verified final Mantel/Spearman block from original script.')

# 11. MANUSCRIPT: TAX4FUN / PAST SIMPER (FIGURE SLOT G) ---------------------------------------------------
# Tax4Fun was run separately using the original SILVA-compatible input and
# its original reference data. No direct ASV-to-SILVA-ID mapping is invented here.
# Expected manuscript input tables (generated externally):
#   input/functional_predition_all(2).xlsx
#   input/pathway_predition_all(2).xlsx
# PAST SIMPER was performed outside R; import the actual PAST result table.
# Do not rerun SIMPER with vegan::simper as a surrogate.
message('ACTION REQUIRED: verify Tax4Fun input/reference and import original PAST SIMPER output.')

# 12. REPRODUCIBILITY ---------------------------------------------------------
writeLines(capture.output(sessionInfo()),file.path(output_dir,'sessionInfo.txt'))
# Software citations: use the exact versions installed in the run environment.
citation_packages <- c('dada2','phyloseq','msa','phangorn','ape','vegan',
                       'ggplot2','dplyr','tidyr','tibble','patchwork','ggrepel')
con <- file(file.path(output_dir,'software_citations.txt'),open='wt')
for (pkg in citation_packages) {
  cat('\n===== ',pkg,' (',as.character(utils::packageVersion(pkg)),') =====\n',
      file=con,sep='')
  cat(paste(capture.output(print(utils::citation(pkg))),collapse='\n'),
      '\n',file=con)
}
cat('\n===== R =====\n',file=con)
cat(paste(capture.output(print(utils::citation())),collapse='\n'),file=con)
close(con)
message('Finished available sections; manuscript-specific blocks flagged as PENDING remain unverified.')


# ============================================================================
# APPENDIX: VERBATIM ORIGINAL FIGURE CODE (COMMENTED, NOT EXECUTED)
# These blocks include interactive commands and historical alternatives.
# They are kept intact as comments, to avoid silently changing the published
# computations. A code audit is required before promoting a block to runnable.
# ============================================================================

# ===== ORIGINAL SOURCE: 01_Fig1_ambientale_CHNS_enzimi.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 11434-13212; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Fig. 1 environmental composite; requires env_full
#| 
#| #   FINAL FIGURE SCRIPT                                                       #
#| #                                                                             #
#| #   Required input object: env_full                                           #
#| #                                                                             #
#| #   Columns:                                                                  #
#| #       E-WI | WI | Control                                                   #
#| #                                                                             #
#| #   Sections:                                                                 #
#| #       A. Physicochemical conditions                                         #
#| #          pH + ORP                                                           #
#| #                                                                             #
#| #       B. Sediment elemental composition                                     #
#| #          C, H, N, S                                                         #
#| #                                                                             #
#| #       C. Extracellular enzymatic activities                                 #
#| #          Leu-MCA + MUF-glu                                                  #
#| #                                                                             #
#| #   T0 = common initial sediment                                              #
#| #                                                                             #
#| #   Sampling positions:                                                       #
#| #       H = H-side sampling position                                          #
#| #       O = O-side sampling position                                          #
#| #                                                                             #
#| #   In E-WI:                                                                  #
#| #       H = cathodic / H2-production region                                   #
#| #       O = anodic / O2-production region                                     #
#| #                                                                             #
#| ###############################################################################
#| 
#| 
#| # =============================================================================
#| # 1. PACKAGES
#| # =============================================================================
#| 
#| library(tidyverse)
#| library(patchwork)
#| 
#| 
#| # =============================================================================
#| # 2. CHECK INPUT DATA
#| # =============================================================================
#| 
#| required_columns <- c(
#|   "Sample",
#|   "Treatment",
#|   "Compartment",
#|   "Time",
#|   
#|   "pH",
#|   "pH_SD",
#|   "ORP",
#|   "ORP_SD",
#|   
#|   "C_mean",
#|   "C_SD",
#|   "H_mean",
#|   "H_SD",
#|   "N_mean",
#|   "N_SD",
#|   "S_mean",
#|   "S_SD",
#|   
#|   "MCA_mean",
#|   "MCA_SD",
#|   "MUF_mean",
#|   "MUF_SD"
#| )
#| 
#| 
#| missing_columns <- setdiff(
#|   required_columns,
#|   colnames(env_full)
#| )
#| 
#| 
#| if(length(missing_columns) > 0) {
#|   
#|   stop(
#|     paste(
#|       "Missing columns in env_full:",
#|       paste(missing_columns, collapse = ", ")
#|     )
#|   )
#| }
#| 
#| 
#| # Check T0
#| 
#| if(sum(env_full$Sample == "T0") != 1) {
#|   
#|   stop(
#|     "env_full must contain exactly one row corresponding to Sample == 'T0'."
#|   )
#| }
#| 
#| 
#| # =============================================================================
#| # 3. PREPARE DATA
#| # =============================================================================
#| 
#| df <- env_full %>%
#|   
#|   mutate(
#|     
#|     Sample = as.character(Sample),
#|     
#|     Treatment = factor(
#|       Treatment,
#|       levels = c(
#|         "E-WI",
#|         "WI",
#|         "Control",
#|         "Initial"
#|       )
#|     ),
#|     
#|     Compartment = factor(
#|       Compartment,
#|       levels = c(
#|         "H",
#|         "O",
#|         "C",
#|         "Initial"
#|       )
#|     ),
#|     
#|     Time = factor(
#|       as.character(Time),
#|       levels = c(
#|         "T0",
#|         "T1",
#|         "T2"
#|       )
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 4. GLOBAL SETTINGS
#| # =============================================================================
#| 
#| treatment_order <- c(
#|   "E-WI",
#|   "WI",
#|   "Control"
#| )
#| 
#| 
#| time_order <- c(
#|   "T0",
#|   "T1",
#|   "T2"
#| )
#| 
#| 
#| # =============================================================================
#| # 5. SAMPLING-POSITION GRAPHICAL CODING
#| #
#| # H = solid line
#| # O = dashed line
#| # C = solid line
#| #
#| # Shape is also used in pH/ORP to reinforce H/O distinction.
#| # =============================================================================
#| 
#| position_linetypes <- c(
#|   "H" = "solid",
#|   "O" = "dashed",
#|   "C" = "solid"
#| )
#| 
#| 
#| position_shapes <- c(
#|   "H" = 16,       # filled circle
#|   "O" = 17,       # filled triangle
#|   "C" = 16
#| )
#| 
#| 
#| # =============================================================================
#| # 6. FIXED ELEMENT COLOURS AND SHAPES
#| #
#| # These are explicitly defined so that plots and legend ALWAYS match.
#| # =============================================================================
#| 
#| element_colors <- c(
#|   "C" = "#F8766D",
#|   "H" = "#7CAE00",
#|   "N" = "#00BFC4",
#|   "S" = "#C77CFF"
#| )
#| 
#| 
#| element_shapes <- c(
#|   "C" = 16,       # circle
#|   "H" = 17,       # triangle
#|   "N" = 15,       # square
#|   "S" = 3         # cross
#| )
#| 
#| 
#| # =============================================================================
#| # 7. FIXED ENZYME COLOURS AND SHAPES
#| # =============================================================================
#| 
#| enzyme_colors <- c(
#|   "Leu-MCA" = "#E76F51",
#|   "MUF-glu" = "#2A9D8F"
#| )
#| 
#| 
#| enzyme_shapes <- c(
#|   "Leu-MCA" = 16,
#|   "MUF-glu" = 17
#| )
#| 
#| 
#| # =============================================================================
#| # 8. COMMON THEME
#| # =============================================================================
#| 
#| theme_environment <- theme_classic(
#|   base_size = 9
#| ) +
#|   
#|   theme(
#|     
#|     axis.text = element_text(
#|       colour = "black",
#|       size = 8
#|     ),
#|     
#|     axis.title.y = element_text(
#|       colour = "black",
#|       size = 8.5
#|     ),
#|     
#|     axis.title.x = element_blank(),
#|     
#|     axis.line = element_line(
#|       linewidth = 0.35
#|     ),
#|     
#|     axis.ticks = element_line(
#|       linewidth = 0.35
#|     ),
#|     
#|     legend.position = "none",
#|     
#|     plot.margin = margin(
#|       3,
#|       5,
#|       3,
#|       5
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 9. FUNCTION TO ADD COMMON T0
#| #
#| # T0 is measured once.
#| #
#| # For graphical representation only:
#| #
#| # E-WI: T0 is duplicated for H and O trajectories
#| # WI:   T0 is duplicated for H and O trajectories
#| # Control: one C trajectory
#| #
#| # No new experimental observation is created.
#| # =============================================================================
#| 
#| add_T0 <- function(
#|     data,
#|     mean_var,
#|     sd_var,
#|     treatment
#| ) {
#|   
#|   # ---------------------------------------------------------------------------
#|   # Existing T1/T2 data
#|   # ---------------------------------------------------------------------------
#|   
#|   d <- data %>%
#|     
#|     filter(
#|       as.character(Treatment) == treatment
#|     ) %>%
#|     
#|     transmute(
#|       Treatment = treatment,
#|       Compartment = as.character(Compartment),
#|       Time = as.character(Time),
#|       Mean = .data[[mean_var]],
#|       SD = .data[[sd_var]]
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Common T0
#|   # ---------------------------------------------------------------------------
#|   
#|   t0 <- data %>%
#|     filter(
#|       Sample == "T0"
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Define trajectories
#|   # ---------------------------------------------------------------------------
#|   
#|   if(treatment == "Control") {
#|     
#|     compartments <- "C"
#|     
#|   } else {
#|     
#|     compartments <- c(
#|       "H",
#|       "O"
#|     )
#|   }
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Duplicate T0 for graphical trajectories
#|   # ---------------------------------------------------------------------------
#|   
#|   t0_rows <- tibble(
#|     
#|     Treatment = treatment,
#|     
#|     Compartment = compartments,
#|     
#|     Time = "T0",
#|     
#|     Mean = rep(
#|       t0[[mean_var]],
#|       length(compartments)
#|     ),
#|     
#|     SD = rep(
#|       t0[[sd_var]],
#|       length(compartments)
#|     )
#|   )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Final object
#|   # ---------------------------------------------------------------------------
#|   
#|   bind_rows(
#|     t0_rows,
#|     d
#|   ) %>%
#|     
#|     mutate(
#|       
#|       Treatment = factor(
#|         Treatment,
#|         levels = treatment_order
#|       ),
#|       
#|       Compartment = factor(
#|         Compartment,
#|         levels = c(
#|           "H",
#|           "O",
#|           "C"
#|         )
#|       ),
#|       
#|       Time = factor(
#|         Time,
#|         levels = time_order
#|       )
#|     ) %>%
#|     
#|     arrange(
#|       Compartment,
#|       Time
#|     )
#| }
#| 
#| 
#| # =============================================================================
#| # 10. pH / ORP PLOT FUNCTION
#| # =============================================================================
#| 
#| make_phys_plot <- function(
#|     mean_var,
#|     sd_var,
#|     treatment,
#|     ylab,
#|     show_x = TRUE,
#|     zero_line = FALSE
#| ) {
#|   
#|   d <- add_T0(
#|     data = df,
#|     mean_var = mean_var,
#|     sd_var = sd_var,
#|     treatment = treatment
#|   )
#|   
#|   
#|   p <- ggplot(
#|     d,
#|     aes(
#|       x = Time,
#|       y = Mean,
#|       group = Compartment,
#|       linetype = Compartment,
#|       shape = Compartment
#|     )
#|   )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # ORP zero reference
#|   # ---------------------------------------------------------------------------
#|   
#|   if(zero_line) {
#|     
#|     p <- p +
#|       
#|       geom_hline(
#|         yintercept = 0,
#|         linetype = "dotted",
#|         linewidth = 0.35,
#|         colour = "grey40"
#|       )
#|   }
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Main plot
#|   # ---------------------------------------------------------------------------
#|   
#|   p <- p +
#|     
#|     geom_line(
#|       linewidth = 0.75,
#|       na.rm = TRUE
#|     ) +
#|     
#|     geom_errorbar(
#|       aes(
#|         ymin = Mean - SD,
#|         ymax = Mean + SD
#|       ),
#|       width = 0.06,
#|       linewidth = 0.35,
#|       na.rm = TRUE
#|     ) +
#|     
#|     geom_point(
#|       size = 2.5,
#|       stroke = 0.7,
#|       na.rm = TRUE
#|     ) +
#|     
#|     scale_linetype_manual(
#|       values = position_linetypes,
#|       drop = FALSE
#|     ) +
#|     
#|     scale_shape_manual(
#|       values = position_shapes,
#|       drop = FALSE
#|     ) +
#|     
#|     labs(
#|       x = NULL,
#|       y = ylab
#|     ) +
#|     
#|     theme_environment
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Remove x labels from upper pH plot
#|   # ---------------------------------------------------------------------------
#|   
#|   if(!show_x) {
#|     
#|     p <- p +
#|       
#|       theme(
#|         axis.text.x = element_blank(),
#|         axis.ticks.x = element_blank()
#|       )
#|   }
#|   
#|   
#|   return(p)
#| }
#| 
#| 
#| # =============================================================================
#| # 11. PHYSICOCHEMICAL BLOCK
#| #
#| # pH and ORP remain separate plots/scales but are visually combined.
#| # =============================================================================
#| 
#| make_phys_block <- function(treatment) {
#|   
#|   p_pH <- make_phys_plot(
#|     mean_var = "pH",
#|     sd_var = "pH_SD",
#|     treatment = treatment,
#|     ylab = "pH",
#|     show_x = FALSE,
#|     zero_line = FALSE
#|   )
#|   
#|   
#|   p_ORP <- make_phys_plot(
#|     mean_var = "ORP",
#|     sd_var = "ORP_SD",
#|     treatment = treatment,
#|     ylab = "ORP (mV)",
#|     show_x = TRUE,
#|     zero_line = TRUE
#|   )
#|   
#|   
#|   p_pH /
#|     p_ORP +
#|     
#|     plot_layout(
#|       heights = c(
#|         1,
#|         1
#|       )
#|     )
#| }
#| 
#| 
#| # =============================================================================
#| # 12. CREATE PHYSICOCHEMICAL PANELS
#| # =============================================================================
#| 
#| phys_EWI <- make_phys_block(
#|   "E-WI"
#| )
#| 
#| phys_WI <- make_phys_block(
#|   "WI"
#| )
#| 
#| phys_Control <- make_phys_block(
#|   "Control"
#| )
#| 
#| 
#| # =============================================================================
#| # 13. CHNS FUNCTION
#| # =============================================================================
#| 
#| make_CHNS <- function(treatment) {
#|   
#|   # ---------------------------------------------------------------------------
#|   # Define trajectories
#|   # ---------------------------------------------------------------------------
#|   
#|   if(treatment == "Control") {
#|     
#|     compartments <- "C"
#|     
#|   } else {
#|     
#|     compartments <- c(
#|       "H",
#|       "O"
#|     )
#|   }
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # T1/T2 observations
#|   # ---------------------------------------------------------------------------
#|   
#|   d <- df %>%
#|     
#|     filter(
#|       as.character(Treatment) == treatment
#|     ) %>%
#|     
#|     transmute(
#|       
#|       Treatment = treatment,
#|       
#|       Compartment = as.character(Compartment),
#|       
#|       Time = as.character(Time),
#|       
#|       C_mean,
#|       C_SD,
#|       
#|       H_mean,
#|       H_SD,
#|       
#|       N_mean,
#|       N_SD,
#|       
#|       S_mean,
#|       S_SD
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Common T0
#|   # ---------------------------------------------------------------------------
#|   
#|   t0 <- df %>%
#|     filter(
#|       Sample == "T0"
#|     )
#|   
#|   
#|   t0_rows <- map_dfr(
#|     
#|     compartments,
#|     
#|     function(comp) {
#|       
#|       tibble(
#|         
#|         Treatment = treatment,
#|         
#|         Compartment = comp,
#|         
#|         Time = "T0",
#|         
#|         C_mean = t0$C_mean,
#|         C_SD = t0$C_SD,
#|         
#|         H_mean = t0$H_mean,
#|         H_SD = t0$H_SD,
#|         
#|         N_mean = t0$N_mean,
#|         N_SD = t0$N_SD,
#|         
#|         S_mean = t0$S_mean,
#|         S_SD = t0$S_SD
#|       )
#|     }
#|   )
#|   
#|   
#|   d <- bind_rows(
#|     t0_rows,
#|     d
#|   )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Convert means to long format
#|   # ---------------------------------------------------------------------------
#|   
#|   means <- d %>%
#|     
#|     select(
#|       Treatment,
#|       Compartment,
#|       Time,
#|       C_mean,
#|       H_mean,
#|       N_mean,
#|       S_mean
#|     ) %>%
#|     
#|     pivot_longer(
#|       cols = c(
#|         C_mean,
#|         H_mean,
#|         N_mean,
#|         S_mean
#|       ),
#|       names_to = "Element",
#|       values_to = "Mean"
#|     ) %>%
#|     
#|     mutate(
#|       Element = str_remove(
#|         Element,
#|         "_mean"
#|       )
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Convert SD to long format
#|   # ---------------------------------------------------------------------------
#|   
#|   sds <- d %>%
#|     
#|     select(
#|       Treatment,
#|       Compartment,
#|       Time,
#|       C_SD,
#|       H_SD,
#|       N_SD,
#|       S_SD
#|     ) %>%
#|     
#|     pivot_longer(
#|       cols = c(
#|         C_SD,
#|         H_SD,
#|         N_SD,
#|         S_SD
#|       ),
#|       names_to = "Element",
#|       values_to = "SD"
#|     ) %>%
#|     
#|     mutate(
#|       Element = str_remove(
#|         Element,
#|         "_SD"
#|       )
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Merge means and SD
#|   # ---------------------------------------------------------------------------
#|   
#|   long <- left_join(
#|     
#|     means,
#|     sds,
#|     
#|     by = c(
#|       "Treatment",
#|       "Compartment",
#|       "Time",
#|       "Element"
#|     )
#|     
#|   ) %>%
#|     
#|     mutate(
#|       
#|       Time = factor(
#|         Time,
#|         levels = time_order
#|       ),
#|       
#|       Element = factor(
#|         Element,
#|         levels = c(
#|           "C",
#|           "H",
#|           "N",
#|           "S"
#|         )
#|       ),
#|       
#|       Compartment = factor(
#|         Compartment,
#|         levels = c(
#|           "H",
#|           "O",
#|           "C"
#|         )
#|       ),
#|       
#|       Series = interaction(
#|         Element,
#|         Compartment,
#|         drop = TRUE
#|       )
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Plot
#|   # ---------------------------------------------------------------------------
#|   
#|   ggplot(
#|     long,
#|     aes(
#|       x = Time,
#|       y = Mean,
#|       group = Series,
#|       colour = Element,
#|       shape = Element,
#|       linetype = Compartment
#|     )
#|   ) +
#|     
#|     geom_line(
#|       linewidth = 0.75,
#|       na.rm = TRUE
#|     ) +
#|     
#|     geom_errorbar(
#|       aes(
#|         ymin = Mean - SD,
#|         ymax = Mean + SD
#|       ),
#|       width = 0.06,
#|       linewidth = 0.30,
#|       na.rm = TRUE
#|     ) +
#|     
#|     geom_point(
#|       size = 2.4,
#|       stroke = 0.7,
#|       na.rm = TRUE
#|     ) +
#|     
#|     scale_colour_manual(
#|       values = element_colors,
#|       breaks = c(
#|         "C",
#|         "H",
#|         "N",
#|         "S"
#|       ),
#|       drop = FALSE
#|     ) +
#|     
#|     scale_shape_manual(
#|       values = element_shapes,
#|       breaks = c(
#|         "C",
#|         "H",
#|         "N",
#|         "S"
#|       ),
#|       drop = FALSE
#|     ) +
#|     
#|     scale_linetype_manual(
#|       values = position_linetypes,
#|       drop = FALSE
#|     ) +
#|     
#|     labs(
#|       x = NULL,
#|       y = "Element content (% DW)"
#|     ) +
#|     
#|     theme_environment
#| }
#| 
#| 
#| # =============================================================================
#| # 14. CREATE CHNS PANELS
#| # =============================================================================
#| 
#| CHNS_EWI <- make_CHNS(
#|   "E-WI"
#| )
#| 
#| CHNS_WI <- make_CHNS(
#|   "WI"
#| )
#| 
#| CHNS_Control <- make_CHNS(
#|   "Control"
#| )
#| 
#| 
#| # =============================================================================
#| # 15. ENZYMATIC ACTIVITY FUNCTION
#| # =============================================================================
#| 
#| make_enzymes <- function(treatment) {
#|   
#|   # ---------------------------------------------------------------------------
#|   # Define trajectories
#|   # ---------------------------------------------------------------------------
#|   
#|   if(treatment == "Control") {
#|     
#|     compartments <- "C"
#|     
#|   } else {
#|     
#|     compartments <- c(
#|       "H",
#|       "O"
#|     )
#|   }
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # T1/T2 observations
#|   # ---------------------------------------------------------------------------
#|   
#|   d <- df %>%
#|     
#|     filter(
#|       as.character(Treatment) == treatment
#|     ) %>%
#|     
#|     transmute(
#|       
#|       Treatment = treatment,
#|       
#|       Compartment = as.character(Compartment),
#|       
#|       Time = as.character(Time),
#|       
#|       MCA_mean,
#|       MCA_SD,
#|       
#|       MUF_mean,
#|       MUF_SD
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Common T0
#|   # ---------------------------------------------------------------------------
#|   
#|   t0 <- df %>%
#|     filter(
#|       Sample == "T0"
#|     )
#|   
#|   
#|   t0_rows <- map_dfr(
#|     
#|     compartments,
#|     
#|     function(comp) {
#|       
#|       tibble(
#|         
#|         Treatment = treatment,
#|         
#|         Compartment = comp,
#|         
#|         Time = "T0",
#|         
#|         MCA_mean = t0$MCA_mean,
#|         MCA_SD = t0$MCA_SD,
#|         
#|         MUF_mean = t0$MUF_mean,
#|         MUF_SD = t0$MUF_SD
#|       )
#|     }
#|   )
#|   
#|   
#|   d <- bind_rows(
#|     t0_rows,
#|     d
#|   )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Means
#|   # ---------------------------------------------------------------------------
#|   
#|   means <- d %>%
#|     
#|     select(
#|       Treatment,
#|       Compartment,
#|       Time,
#|       MCA_mean,
#|       MUF_mean
#|     ) %>%
#|     
#|     pivot_longer(
#|       cols = c(
#|         MCA_mean,
#|         MUF_mean
#|       ),
#|       names_to = "Enzyme",
#|       values_to = "Mean"
#|     ) %>%
#|     
#|     mutate(
#|       
#|       Enzyme = recode(
#|         Enzyme,
#|         "MCA_mean" = "Leu-MCA",
#|         "MUF_mean" = "MUF-glu"
#|       )
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Standard deviations
#|   # ---------------------------------------------------------------------------
#|   
#|   sds <- d %>%
#|     
#|     select(
#|       Treatment,
#|       Compartment,
#|       Time,
#|       MCA_SD,
#|       MUF_SD
#|     ) %>%
#|     
#|     pivot_longer(
#|       cols = c(
#|         MCA_SD,
#|         MUF_SD
#|       ),
#|       names_to = "Enzyme",
#|       values_to = "SD"
#|     ) %>%
#|     
#|     mutate(
#|       
#|       Enzyme = recode(
#|         Enzyme,
#|         "MCA_SD" = "Leu-MCA",
#|         "MUF_SD" = "MUF-glu"
#|       )
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Merge
#|   # ---------------------------------------------------------------------------
#|   
#|   long <- left_join(
#|     
#|     means,
#|     sds,
#|     
#|     by = c(
#|       "Treatment",
#|       "Compartment",
#|       "Time",
#|       "Enzyme"
#|     )
#|     
#|   ) %>%
#|     
#|     mutate(
#|       
#|       Time = factor(
#|         Time,
#|         levels = time_order
#|       ),
#|       
#|       Enzyme = factor(
#|         Enzyme,
#|         levels = c(
#|           "Leu-MCA",
#|           "MUF-glu"
#|         )
#|       ),
#|       
#|       Compartment = factor(
#|         Compartment,
#|         levels = c(
#|           "H",
#|           "O",
#|           "C"
#|         )
#|       ),
#|       
#|       Series = interaction(
#|         Enzyme,
#|         Compartment,
#|         drop = TRUE
#|       )
#|     )
#|   
#|   
#|   # ---------------------------------------------------------------------------
#|   # Plot
#|   # ---------------------------------------------------------------------------
#|   
#|   ggplot(
#|     long,
#|     aes(
#|       x = Time,
#|       y = Mean,
#|       group = Series,
#|       colour = Enzyme,
#|       shape = Enzyme,
#|       linetype = Compartment
#|     )
#|   ) +
#|     
#|     geom_hline(
#|       yintercept = 0,
#|       linetype = "dotted",
#|       linewidth = 0.30,
#|       colour = "grey40"
#|     ) +
#|     
#|     geom_line(
#|       linewidth = 0.75,
#|       na.rm = TRUE
#|     ) +
#|     
#|     geom_errorbar(
#|       aes(
#|         ymin = Mean - SD,
#|         ymax = Mean + SD
#|       ),
#|       width = 0.06,
#|       linewidth = 0.30,
#|       na.rm = TRUE
#|     ) +
#|     
#|     geom_point(
#|       size = 2.4,
#|       stroke = 0.7,
#|       na.rm = TRUE
#|     ) +
#|     
#|     scale_colour_manual(
#|       values = enzyme_colors,
#|       breaks = c(
#|         "Leu-MCA",
#|         "MUF-glu"
#|       ),
#|       drop = FALSE
#|     ) +
#|     
#|     scale_shape_manual(
#|       values = enzyme_shapes,
#|       breaks = c(
#|         "Leu-MCA",
#|         "MUF-glu"
#|       ),
#|       drop = FALSE
#|     ) +
#|     
#|     scale_linetype_manual(
#|       values = position_linetypes,
#|       drop = FALSE
#|     ) +
#|     
#|     labs(
#|       x = NULL,
#|       
#|       y = expression(
#|         "Activity (nmol g"^{-1}*" h"^{-1}*")"
#|       )
#|     ) +
#|     
#|     theme_environment
#| }
#| 
#| 
#| # =============================================================================
#| # 16. CREATE ENZYME PANELS
#| # =============================================================================
#| 
#| ENZ_EWI <- make_enzymes(
#|   "E-WI"
#| )
#| 
#| ENZ_WI <- make_enzymes(
#|   "WI"
#| )
#| 
#| ENZ_Control <- make_enzymes(
#|   "Control"
#| )
#| 
#| 
#| # =============================================================================
#| # 17. COLUMN HEADERS
#| # =============================================================================
#| 
#| make_header <- function(label) {
#|   
#|   ggplot() +
#|     
#|     annotate(
#|       "text",
#|       x = 0.5,
#|       y = 0.5,
#|       label = label,
#|       fontface = "bold",
#|       size = 4.3
#|     ) +
#|     
#|     xlim(
#|       0,
#|       1
#|     ) +
#|     
#|     ylim(
#|       0,
#|       1
#|     ) +
#|     
#|     theme_void()
#| }
#| 
#| 
#| header_EWI <- make_header(
#|   "E-WI"
#| )
#| 
#| header_WI <- make_header(
#|   "WI"
#| )
#| 
#| header_Control <- make_header(
#|   "Control"
#| )
#| 
#| 
#| # =============================================================================
#| # 18. SECTION TITLES
#| # =============================================================================
#| 
#| make_section_title <- function(label) {
#|   
#|   ggplot() +
#|     
#|     annotate(
#|       "text",
#|       x = 0.01,
#|       y = 0.5,
#|       label = label,
#|       hjust = 0,
#|       vjust = 0.5,
#|       fontface = "bold",
#|       size = 3.7
#|     ) +
#|     
#|     xlim(
#|       0,
#|       1
#|     ) +
#|     
#|     ylim(
#|       0,
#|       1
#|     ) +
#|     
#|     theme_void() +
#|     
#|     theme(
#|       plot.margin = margin(
#|         0,
#|         0,
#|         0,
#|         4
#|       )
#|     )
#| }
#| 
#| 
#| title_A <- make_section_title(
#|   "A. Physicochemical conditions"
#| )
#| 
#| title_B <- make_section_title(
#|   "B. Sediment elemental composition"
#| )
#| 
#| title_C <- make_section_title(
#|   "C. Extracellular enzymatic activities"
#| )
#| 
#| 
#| # =============================================================================
#| # 19. BUILD MAIN DATA ROWS
#| # =============================================================================
#| 
#| row_A <- (
#|   phys_EWI |
#|     phys_WI |
#|     phys_Control
#| )
#| 
#| 
#| row_B <- (
#|   CHNS_EWI |
#|     CHNS_WI |
#|     CHNS_Control
#| )
#| 
#| 
#| row_C <- (
#|   ENZ_EWI |
#|     ENZ_WI |
#|     ENZ_Control
#| )
#| 
#| 
#| # =============================================================================
#| # 20. COMMON LEGEND
#| #
#| # IMPORTANT:
#| #
#| # This legend does NOT use automatic ggplot colour scales.
#| #
#| # Colours and shapes are explicitly taken from:
#| #
#| #   element_colors
#| #   element_shapes
#| #   enzyme_colors
#| #   enzyme_shapes
#| #
#| # Therefore the legend and the plots cannot become mismatched.
#| # =============================================================================
#| 
#| 
#| # -----------------------------------------------------------------------------
#| # Element legend data
#| # -----------------------------------------------------------------------------
#| 
#| legend_elements <- tibble(
#|   
#|   Element = factor(
#|     c(
#|       "C",
#|       "H",
#|       "N",
#|       "S"
#|     ),
#|     levels = c(
#|       "C",
#|       "H",
#|       "N",
#|       "S"
#|     )
#|   ),
#|   
#|   x = c(
#|     1.30,
#|     2.05,
#|     2.80,
#|     3.55
#|   ),
#|   
#|   y = 2
#| )
#| 
#| 
#| # -----------------------------------------------------------------------------
#| # Enzyme legend data
#| # -----------------------------------------------------------------------------
#| 
#| legend_enzymes <- tibble(
#|   
#|   Enzyme = factor(
#|     c(
#|       "Leu-MCA",
#|       "MUF-glu"
#|     ),
#|     levels = c(
#|       "Leu-MCA",
#|       "MUF-glu"
#|     )
#|   ),
#|   
#|   x = c(
#|     1.80,
#|     3.20
#|   ),
#|   
#|   y = 1
#| )
#| 
#| 
#| # =============================================================================
#| # 21. DRAW COMMON LEGEND
#| # =============================================================================
#| 
#| legend_plot <- ggplot() +
#|   
#|   
#|   # ---------------------------------------------------------------------------
#| # SAMPLING POSITION
#| # ---------------------------------------------------------------------------
#| 
#| annotate(
#|   "text",
#|   x = 0,
#|   y = 3,
#|   label = "Sampling position:",
#|   hjust = 0,
#|   fontface = "bold",
#|   size = 3.2
#| ) +
#|   
#|   
#|   # H
#|   
#|   annotate(
#|     "segment",
#|     x = 1.40,
#|     xend = 1.90,
#|     y = 3,
#|     yend = 3,
#|     linewidth = 0.75,
#|     linetype = "solid"
#|   ) +
#|   
#|   annotate(
#|     "point",
#|     x = 1.65,
#|     y = 3,
#|     shape = 16,
#|     size = 2.4
#|   ) +
#|   
#|   annotate(
#|     "text",
#|     x = 2.00,
#|     y = 3,
#|     label = "H",
#|     hjust = 0,
#|     size = 3.2
#|   ) +
#|   
#|   
#|   # O
#|   
#|   annotate(
#|     "segment",
#|     x = 2.50,
#|     xend = 3.00,
#|     y = 3,
#|     yend = 3,
#|     linewidth = 0.75,
#|     linetype = "dashed"
#|   ) +
#|   
#|   annotate(
#|     "point",
#|     x = 2.75,
#|     y = 3,
#|     shape = 17,
#|     size = 2.4
#|   ) +
#|   
#|   annotate(
#|     "text",
#|     x = 3.10,
#|     y = 3,
#|     label = "O",
#|     hjust = 0,
#|     size = 3.2
#|   ) +
#|   
#|   
#|   # ---------------------------------------------------------------------------
#| # ELEMENTS
#| # ---------------------------------------------------------------------------
#| 
#| annotate(
#|   "text",
#|   x = 0,
#|   y = 2,
#|   label = "Element:",
#|   hjust = 0,
#|   fontface = "bold",
#|   size = 3.2
#| ) +
#|   
#|   
#|   # C
#|   
#|   annotate(
#|     "point",
#|     x = 1.30,
#|     y = 2,
#|     colour = element_colors["C"],
#|     shape = element_shapes["C"],
#|     size = 2.5
#|   ) +
#|   
#|   annotate(
#|     "text",
#|     x = 1.48,
#|     y = 2,
#|     label = "C",
#|     hjust = 0,
#|     size = 3.2
#|   ) +
#|   
#|   
#|   # H
#|   
#|   annotate(
#|     "point",
#|     x = 2.05,
#|     y = 2,
#|     colour = element_colors["H"],
#|     shape = element_shapes["H"],
#|     size = 2.5
#|   ) +
#|   
#|   annotate(
#|     "text",
#|     x = 2.23,
#|     y = 2,
#|     label = "H",
#|     hjust = 0,
#|     size = 3.2
#|   ) +
#|   
#|   
#|   # N
#|   
#|   annotate(
#|     "point",
#|     x = 2.80,
#|     y = 2,
#|     colour = element_colors["N"],
#|     shape = element_shapes["N"],
#|     size = 2.5
#|   ) +
#|   
#|   annotate(
#|     "text",
#|     x = 2.98,
#|     y = 2,
#|     label = "N",
#|     hjust = 0,
#|     size = 3.2
#|   ) +
#|   
#|   
#|   # S
#|   
#|   annotate(
#|     "point",
#|     x = 3.55,
#|     y = 2,
#|     colour = element_colors["S"],
#|     shape = element_shapes["S"],
#|     size = 2.5
#|   ) +
#|   
#|   annotate(
#|     "text",
#|     x = 3.73,
#|     y = 2,
#|     label = "S",
#|     hjust = 0,
#|     size = 3.2
#|   ) +
#|   
#|   
#|   # ---------------------------------------------------------------------------
#| # ENZYMATIC ACTIVITY
#| # ---------------------------------------------------------------------------
#| 
#| annotate(
#|   "text",
#|   x = 0,
#|   y = 1,
#|   label = "Enzymatic activity:",
#|   hjust = 0,
#|   fontface = "bold",
#|   size = 3.2
#| ) +
#|   
#|   
#|   # Leu-MCA
#|   
#|   annotate(
#|     "point",
#|     x = 1.80,
#|     y = 1,
#|     colour = enzyme_colors["Leu-MCA"],
#|     shape = enzyme_shapes["Leu-MCA"],
#|     size = 2.5
#|   ) +
#|   
#|   annotate(
#|     "text",
#|     x = 1.98,
#|     y = 1,
#|     label = "Leu-MCA",
#|     hjust = 0,
#|     size = 3.2
#|   ) +
#|   
#|   
#|   # MUF-glu
#|   
#|   annotate(
#|     "point",
#|     x = 3.20,
#|     y = 1,
#|     colour = enzyme_colors["MUF-glu"],
#|     shape = enzyme_shapes["MUF-glu"],
#|     size = 2.5
#|   ) +
#|   
#|   annotate(
#|     "text",
#|     x = 3.38,
#|     y = 1,
#|     label = "MUF-glu",
#|     hjust = 0,
#|     size = 3.2
#|   ) +
#|   
#|   
#|   # ---------------------------------------------------------------------------
#| # Legend layout
#| # ---------------------------------------------------------------------------
#| 
#| coord_cartesian(
#|   xlim = c(
#|     0,
#|     5
#|   ),
#|   ylim = c(
#|     0.5,
#|     3.5
#|   ),
#|   clip = "off"
#| ) +
#|   
#|   theme_void() +
#|   
#|   theme(
#|     
#|     plot.margin = margin(
#|       2,
#|       10,
#|       2,
#|       10
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 22. TREATMENT HEADERS
#| # =============================================================================
#| 
#| headers <- (
#|   header_EWI |
#|     header_WI |
#|     header_Control
#| )
#| 
#| 
#| # =============================================================================
#| # 23. SECTION A
#| # =============================================================================
#| 
#| section_A <- (
#|   
#|   title_A /
#|     
#|     row_A
#|   
#| ) +
#|   
#|   plot_layout(
#|     heights = c(
#|       0.075,
#|       1
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 24. SECTION B
#| # =============================================================================
#| 
#| section_B <- (
#|   
#|   title_B /
#|     
#|     row_B
#|   
#| ) +
#|   
#|   plot_layout(
#|     heights = c(
#|       0.075,
#|       1
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 25. SECTION C
#| # =============================================================================
#| 
#| section_C <- (
#|   
#|   title_C /
#|     
#|     row_C
#|   
#| ) +
#|   
#|   plot_layout(
#|     heights = c(
#|       0.075,
#|       1
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 26. FINAL FIGURE
#| # =============================================================================
#| 
#| Figure_environment_FINAL <- (
#|   
#|   headers /
#|     
#|     section_A /
#|     
#|     section_B /
#|     
#|     section_C /
#|     
#|     legend_plot
#|   
#| ) +
#|   
#|   plot_layout(
#|     
#|     heights = c(
#|       0.08,   # treatment headers
#|       1.55,   # physicochemical block
#|       1.00,   # CHNS
#|       1.00,   # enzymatic activity
#|       0.30    # common legend
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 27. DISPLAY
#| # =============================================================================
#| 
#| print(
#|   Figure_environment_FINAL
#| )
#| 
#| 
#| # =============================================================================
#| # 28. EXPORT
#| # =============================================================================
#| 
#| ggsave(
#|   filename = "Figure_environment_FINAL.pdf",
#|   plot = Figure_environment_FINAL,
#|   width = 13,
#|   height = 10.5,
#|   units = "in"
#| )
#| 
#| 
#| ggsave(
#|   filename = "Figure_environment_FINAL.svg",
#|   plot = Figure_environment_FINAL,
#|   width = 13,
#|   height = 10.5,
#|   units = "in"
#| )
#| 
#| 
#| ggsave(
#|   filename = "Figure_environment_FINAL.png",
#|   plot = Figure_environment_FINAL,
#|   width = 13,
#|   height = 10.5,
#|   units = "in",
#|   dpi = 600
#| )

# ===== ORIGINAL SOURCE: 02_Fig1_ambientale_alternativa.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 9303-9909; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Alternative environmental figure code; archival
#| 
#| # FIGURE - TEMPORAL EVOLUTION OF ENVIRONMENTAL AND ENZYMATIC VARIABLES
#| # Dataset: env_full
#| # T0 = common initial sediment condition
#| ###############################################################################
#| 
#| library(tidyverse)
#| library(patchwork)
#| 
#| # ============================================================
#| # 1. PREPARE DATA
#| # ============================================================
#| 
#| plot_data <- env_full %>%
#|   mutate(Group = case_when(
#|     Treatment == "E-WI" & Compartment == "H" ~ "E-WI H",
#|     Treatment == "E-WI" & Compartment == "O" ~ "E-WI O",
#|     Treatment == "WI"   & Compartment == "H" ~ "WI H",
#|     Treatment == "WI"   & Compartment == "O" ~ "WI O",
#|     Treatment == "Control" ~ "Control",
#|     TRUE ~ NA_character_
#|   ))
#| 
#| group_levels <- c("E-WI H", "E-WI O", "WI H", "WI O", "Control")
#| 
#| # ============================================================
#| # 2. FUNCTION TO BUILD EACH TEMPORAL PANEL
#| # ============================================================
#| 
#| make_panel <- function(mean_var, sd_var, title, ylab, zero_line = FALSE) {
#|   
#|   # T1-T2 observations
#|   d <- plot_data %>%
#|     filter(!is.na(Group)) %>%
#|     transmute(Group, Time, Mean = .data[[mean_var]], SD = .data[[sd_var]])
#|   
#|   # Common T0 duplicated only for graphical representation
#|   t0 <- env_full %>%
#|     filter(Sample == "T0") %>%
#|     transmute(Mean = .data[[mean_var]], SD = .data[[sd_var]])
#|   
#|   t0_plot <- tibble(
#|     Group = group_levels,
#|     Time = factor("T0", levels = c("T0", "T1", "T2")),
#|     Mean = rep(t0$Mean, length(group_levels)),
#|     SD = rep(t0$SD, length(group_levels))
#|   )
#|   
#|   d <- bind_rows(t0_plot, d) %>%
#|     mutate(
#|       Group = factor(Group, levels = group_levels),
#|       Time = factor(as.character(Time), levels = c("T0", "T1", "T2"))
#|     )
#|   
#|   p <- ggplot(d, aes(x = Time, y = Mean, group = Group,
#|                      colour = Group, shape = Group)) +
#|     
#|     geom_line(linewidth = 0.75, na.rm = TRUE) +
#|     
#|     geom_errorbar(
#|       aes(ymin = Mean - SD, ymax = Mean + SD),
#|       width = 0.08,
#|       linewidth = 0.45,
#|       na.rm = TRUE
#|     ) +
#|     
#|     geom_point(size = 2.8, stroke = 0.8, na.rm = TRUE) +
#|     
#|     labs(title = title, x = NULL, y = ylab, colour = NULL, shape = NULL) +
#|     
#|     theme_classic(base_size = 10) +
#|     
#|     theme(
#|       plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
#|       axis.text = element_text(colour = "black"),
#|       axis.title.y = element_text(size = 9),
#|       legend.position = "none"
#|     )
#|   
#|   if(zero_line)
#|     p <- p + geom_hline(yintercept = 0, linetype = "dashed",
#|                         linewidth = 0.35, inherit.aes = FALSE)
#|   
#|   return(p)
#| }
#| 
#| # ============================================================
#| # 3. PHYSICOCHEMICAL VARIABLES
#| # ============================================================
#| 
#| p_pH <- make_panel(
#|   "pH", "pH_SD",
#|   "A. pH",
#|   "pH"
#| )
#| 
#| p_ORP <- make_panel(
#|   "ORP", "ORP_SD",
#|   "B. Redox potential",
#|   "ORP (mV)",
#|   zero_line = TRUE
#| )
#| 
#| # ============================================================
#| # 4. ELEMENTAL COMPOSITION
#| # ============================================================
#| 
#| p_C <- make_panel(
#|   "C_mean", "C_SD",
#|   "C. Carbon",
#|   "C (% DW)"
#| )
#| 
#| p_H <- make_panel(
#|   "H_mean", "H_SD",
#|   "D. Hydrogen",
#|   "H (% DW)"
#| )
#| 
#| p_N <- make_panel(
#|   "N_mean", "N_SD",
#|   "E. Nitrogen",
#|   "N (% DW)"
#| )
#| 
#| p_S <- make_panel(
#|   "S_mean", "S_SD",
#|   "F. Sulfur",
#|   "S (% DW)"
#| )
#| 
#| # ============================================================
#| # 5. EXTRACELLULAR ENZYMATIC ACTIVITIES
#| # ============================================================
#| 
#| p_MCA <- make_panel(
#|   "MCA_mean", "MCA_SD",
#|   "G. Leucine aminopeptidase",
#|   expression("Leu-MCA (nmol g"^{-1}*" h"^{-1}*")")
#| )
#| 
#| p_MUF <- make_panel(
#|   "MUF_mean", "MUF_SD",
#|   "H. β-D-glucosidase",
#|   expression("MUF-glu (nmol g"^{-1}*" h"^{-1}*")")
#| )
#| 
#| # ============================================================
#| # 6. ADD ONE COMMON LEGEND
#| # ============================================================
#| 
#| # Keep legend only in first panel
#| p_pH <- p_pH +
#|   theme(legend.position = "bottom") +
#|   guides(
#|     colour = guide_legend(nrow = 1),
#|     shape  = guide_legend(nrow = 1)
#|   )
#| 
#| # ============================================================
#| # 7. COMBINE PANELS
#| # ============================================================
#| 
#| Figure_environment <- (
#|   p_pH | p_ORP | p_C | p_H
#| ) / (
#|   p_N | p_S | p_MCA | p_MUF
#| ) +
#|   plot_layout(guides = "collect") &
#|   theme(legend.position = "bottom")
#| 
#| print(Figure_environment)
#| 
#| # ============================================================
#| # 8. EXPORT
#| # ============================================================
#| 
#| ggsave("Figure_environment_full.pdf", Figure_environment,
#|        width = 15, height = 8, units = "in")
#| 
#| ggsave("Figure_environment_full.svg", Figure_environment,
#|        width = 15, height = 8, units = "in")
#| 
#| ggsave("Figure_environment_full.png", Figure_environment,
#|        width = 15, height = 8, units = "in", dpi = 600)
#| ###############################################################################
#| ###############################################################################
#| ###############################################################################
#| # FIGURE 1B-D
#| # TEMPORAL EVOLUTION OF PHYSICOCHEMICAL CONDITIONS,
#| # SEDIMENT ELEMENTAL COMPOSITION AND ENZYMATIC ACTIVITIES
#| #
#| # Input object: env_full
#| #
#| # T0 = common initial sediment condition
#| # V1 = E-WI
#| # V2 = WI
#| # V3 = Control
#| #
#| # H = cathodic / H2-production region
#| # O = anodic / O2-production region
#| ###############################################################################
#| 
#| library(tidyverse)
#| library(patchwork)
#| 
#| 
#| # =============================================================================
#| # 1. PREPARE DATA
#| # =============================================================================
#| 
#| plot_data <- env_full %>%
#|   mutate(
#|     Group = case_when(
#|       Treatment == "E-WI" & Compartment == "H" ~ "E-WI H",
#|       Treatment == "E-WI" & Compartment == "O" ~ "E-WI O",
#|       Treatment == "WI"   & Compartment == "H" ~ "WI H",
#|       Treatment == "WI"   & Compartment == "O" ~ "WI O",
#|       Treatment == "Control"                   ~ "Control",
#|       TRUE                                     ~ NA_character_
#|     )
#|   )
#| 
#| group_levels <- c(
#|   "E-WI H",
#|   "E-WI O",
#|   "WI H",
#|   "WI O",
#|   "Control"
#| )
#| 
#| 
#| # =============================================================================
#| # 2. GENERAL FUNCTION
#| # =============================================================================
#| 
#| make_panel <- function(mean_var,
#|                        sd_var,
#|                        title,
#|                        ylab,
#|                        zero_line = FALSE) {
#|   
#|   # -----------------------------------------------------------
#|   # T1 and T2
#|   # -----------------------------------------------------------
#|   
#|   d <- plot_data %>%
#|     filter(!is.na(Group)) %>%
#|     transmute(
#|       Group,
#|       Time,
#|       Mean = .data[[mean_var]],
#|       SD   = .data[[sd_var]]
#|     )
#|   
#|   
#|   # -----------------------------------------------------------
#|   # Common T0
#|   #
#|   # T0 is duplicated here ONLY for graphical representation.
#|   # It remains a single common experimental baseline in env_full.
#|   # -----------------------------------------------------------
#|   
#|   t0 <- env_full %>%
#|     filter(Sample == "T0") %>%
#|     transmute(
#|       Mean = .data[[mean_var]],
#|       SD   = .data[[sd_var]]
#|     )
#|   
#|   t0_plot <- tibble(
#|     Group = group_levels,
#|     
#|     Time = factor(
#|       rep("T0", length(group_levels)),
#|       levels = c("T0", "T1", "T2")
#|     ),
#|     
#|     Mean = rep(t0$Mean, length(group_levels)),
#|     SD   = rep(t0$SD, length(group_levels))
#|   )
#|   
#|   
#|   # -----------------------------------------------------------
#|   # Combine T0 + experimental observations
#|   # -----------------------------------------------------------
#|   
#|   d <- bind_rows(t0_plot, d) %>%
#|     mutate(
#|       Group = factor(
#|         Group,
#|         levels = group_levels
#|       ),
#|       
#|       Time = factor(
#|         as.character(Time),
#|         levels = c("T0", "T1", "T2")
#|       )
#|     )
#|   
#|   
#|   # -----------------------------------------------------------
#|   # Plot
#|   # -----------------------------------------------------------
#|   
#|   p <- ggplot(
#|     d,
#|     aes(
#|       x = Time,
#|       y = Mean,
#|       group = Group,
#|       colour = Group,
#|       shape = Group
#|     )
#|   ) +
#|     
#|     geom_line(
#|       linewidth = 0.75,
#|       na.rm = TRUE
#|     ) +
#|     
#|     geom_errorbar(
#|       aes(
#|         ymin = Mean - SD,
#|         ymax = Mean + SD
#|       ),
#|       width = 0.08,
#|       linewidth = 0.45,
#|       na.rm = TRUE
#|     ) +
#|     
#|     geom_point(
#|       size = 2.8,
#|       stroke = 0.8,
#|       na.rm = TRUE
#|     ) +
#|     
#|     labs(
#|       title = title,
#|       x = NULL,
#|       y = ylab,
#|       colour = NULL,
#|       shape = NULL
#|     ) +
#|     
#|     theme_classic(
#|       base_size = 10
#|     ) +
#|     
#|     theme(
#|       
#|       plot.title = element_text(
#|         face = "bold",
#|         size = 10,
#|         hjust = 0.5
#|       ),
#|       
#|       axis.text = element_text(
#|         colour = "black"
#|       ),
#|       
#|       axis.title.y = element_text(
#|         size = 9
#|       ),
#|       
#|       legend.position = "none"
#|     )
#|   
#|   
#|   # -----------------------------------------------------------
#|   # Optional zero reference
#|   # -----------------------------------------------------------
#|   
#|   if(zero_line) {
#|     
#|     p <- p +
#|       geom_hline(
#|         yintercept = 0,
#|         linetype = "dashed",
#|         linewidth = 0.35,
#|         inherit.aes = FALSE
#|       )
#|   }
#|   
#|   return(p)
#| }
#| 
#| 
#| # =============================================================================
#| # 3. PANEL B
#| # PHYSICOCHEMICAL CONDITIONS
#| # =============================================================================
#| 
#| p_pH <- make_panel(
#|   mean_var = "pH",
#|   sd_var   = "pH_SD",
#|   title    = "pH",
#|   ylab     = "pH"
#| )
#| 
#| 
#| p_ORP <- make_panel(
#|   mean_var = "ORP",
#|   sd_var   = "ORP_SD",
#|   title    = "Redox potential",
#|   ylab     = "ORP (mV)",
#|   zero_line = TRUE
#| )
#| 
#| 
#| # =============================================================================
#| # 4. PANEL C
#| # SEDIMENT ELEMENTAL COMPOSITION
#| # =============================================================================
#| 
#| p_C <- make_panel(
#|   mean_var = "C_mean",
#|   sd_var   = "C_SD",
#|   title    = "Carbon",
#|   ylab     = "C (% DW)"
#| )
#| 
#| 
#| p_H <- make_panel(
#|   mean_var = "H_mean",
#|   sd_var   = "H_SD",
#|   title    = "Hydrogen",
#|   ylab     = "H (% DW)"
#| )
#| 
#| 
#| p_N <- make_panel(
#|   mean_var = "N_mean",
#|   sd_var   = "N_SD",
#|   title    = "Nitrogen",
#|   ylab     = "N (% DW)"
#| )
#| 
#| 
#| p_S <- make_panel(
#|   mean_var = "S_mean",
#|   sd_var   = "S_SD",
#|   title    = "Sulfur",
#|   ylab     = "S (% DW)"
#| )
#| 
#| 
#| # =============================================================================
#| # 5. PANEL D
#| # EXTRACELLULAR ENZYMATIC ACTIVITIES
#| # =============================================================================
#| 
#| p_MCA <- make_panel(
#|   mean_var = "MCA_mean",
#|   sd_var   = "MCA_SD",
#|   title    = "Leu-MCA",
#|   ylab     = expression(
#|     "Activity (nmol g"^{-1}*" h"^{-1}*")"
#|   )
#| )
#| 
#| 
#| p_MUF <- make_panel(
#|   mean_var = "MUF_mean",
#|   sd_var   = "MUF_SD",
#|   title    = "MUF-glu",
#|   ylab     = expression(
#|     "Activity (nmol g"^{-1}*" h"^{-1}*")"
#|   )
#| )
#| 
#| 
#| # =============================================================================
#| # 6. BUILD INDIVIDUAL BLOCKS
#| # =============================================================================
#| 
#| # ---------------------------------------------------------------------------
#| # B. Physicochemical conditions
#| # ---------------------------------------------------------------------------
#| 
#| panel_B <- (
#|   p_pH |
#|     p_ORP
#| ) +
#|   plot_annotation(
#|     title = "B. Physicochemical conditions"
#|   )
#| 
#| 
#| # ---------------------------------------------------------------------------
#| # C. Sediment elemental composition
#| # ---------------------------------------------------------------------------
#| 
#| panel_C <- (
#|   p_C |
#|     p_H |
#|     p_N |
#|     p_S
#| ) +
#|   plot_annotation(
#|     title = "C. Sediment elemental composition"
#|   )
#| 
#| 
#| # ---------------------------------------------------------------------------
#| # D. Extracellular enzymatic activities
#| # ---------------------------------------------------------------------------
#| 
#| panel_D <- (
#|   p_MCA |
#|     p_MUF
#| ) +
#|   plot_annotation(
#|     title = "D. Extracellular enzymatic activities"
#|   )
#| 
#| 
#| # =============================================================================
#| # 7. RESTORE ONE LEGEND
#| # =============================================================================
#| 
#| # The legend is activated in one plot.
#| # patchwork will collect it and use it as a common legend.
#| 
#| p_pH <- p_pH +
#|   theme(
#|     legend.position = "bottom"
#|   ) +
#|   guides(
#|     colour = guide_legend(
#|       nrow = 1,
#|       byrow = TRUE
#|     ),
#|     shape = guide_legend(
#|       nrow = 1,
#|       byrow = TRUE
#|     )
#|   )
#| 
#| 
#| # Rebuild panel B after restoring legend
#| 
#| panel_B <- (
#|   p_pH |
#|     p_ORP
#| ) +
#|   plot_annotation(
#|     title = "B. Physicochemical conditions"
#|   )
#| 
#| 
#| # =============================================================================
#| # 8. COMPLETE FIGURE B-D
#| # =============================================================================
#| 
#| Figure1_BD <- (
#|   panel_B /
#|     panel_C /
#|     panel_D
#| ) +
#|   plot_layout(
#|     heights = c(
#|       1.0,   # pH + ORP
#|       1.05,  # CHNS
#|       1.0    # enzymes
#|     ),
#|     guides = "collect"
#|   ) &
#|   theme(
#|     legend.position = "bottom"
#|   )
#| 
#| 
#| # =============================================================================
#| # 9. DISPLAY
#| # =============================================================================
#| 
#| print(Figure1_BD)
#| 
#| 
#| # =============================================================================
#| # 10. EXPORT
#| # =============================================================================
#| 
#| ggsave(
#|   filename = "Figure1_BD_environmental_enzymatic.pdf",
#|   plot = Figure1_BD,
#|   width = 14,
#|   height = 10,
#|   units = "in"
#| )
#| 
#| 
#| ggsave(
#|   filename = "Figure1_BD_environmental_enzymatic.svg",
#|   plot = Figure1_BD,
#|   width = 14,
#|   height = 10,
#|   units = "in"
#| )
#| 
#| 
#| ggsave(
#|   filename = "Figure1_BD_environmental_enzymatic.png",
#|   plot = Figure1_BD,
#|   width = 14,
#|   height = 10,
#|   units = "in",

# ===== ORIGINAL SOURCE: 03_Fig2_alpha_diversity.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 6422-6489; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Alpha diversity (original Fig.1 -> manuscript Fig.2)
#| 
#| # 2. FIGURE 1 — ALPHA DIVERSITY
#| # ============================================================
#| 
#| alpha_div <- estimate_richness(
#|   prok_data,
#|   measures = c(
#|     "Observed",
#|     "Shannon",
#|     "Simpson"
#|   )
#| ) %>%
#|   rownames_to_column("Sample") %>%
#|   left_join( meta %>% select( Sample, treatment, time, electrode ), by = "Sample" )
#| 
#| alpha_div$Sample <- factor(
#|   alpha_div$Sample,
#|   levels = sample_levels
#| )
#| 
#| # Print numerical values
#| cat("\n================ ALPHA DIVERSITY ================\n")
#| 
#| alpha_div %>%
#|   arrange(Sample) %>%
#|   as_tibble() %>%
#|   print(n = Inf)
#| 
#| # Save complete table
#| write.csv( alpha_div, "Table_alpha_diversity.csv", row.names = FALSE )
#| 
#| # Figure: Observed + Shannon
#| # Simpson remains available in supplementary table
#| alpha_long <- alpha_div %>%
#|   select( Sample, treatment, time, electrode, Observed, Shannon ) %>%
#|   pivot_longer( cols = c( Observed, Shannon ), names_to = "Metric", values_to = "Value" )
#| 
#| p_alpha <- ggplot(
#|   alpha_long,
#|   aes( x = Sample, y = Value, colour = treatment, shape = electrode )
#| ) +
#|   geom_point( size = 4 ) +
#|   facet_wrap( ~Metric, scales = "free_y", ncol = 1 ) +
#|   labs( x = NULL, y = NULL, colour = "Treatment", shape = "Electrode" ) +
#|   theme_classic( base_size = 11 ) +
#|   theme(
#|     axis.text.x = element_text(
#|       angle = 45,
#|       hjust = 1,
#|       face = "bold",
#|       colour = "black"
#|     ),
#|     strip.background = element_blank(),
#|     strip.text = element_text(
#|       face = "bold"
#|     ),
#|     legend.position = "right"
#|   )
#| 
#| p_alpha
#| 
#| save_fig(
#|   p_alpha,
#|   "Figure1_Alpha_diversity",
#|   8,
#|   7
#| )
#| 
#| # ============================================================

# ===== ORIGINAL SOURCE: 04_Fig3_UniFrac.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 6490-6660; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Weighted and unweighted UniFrac (original Fig.2 -> manuscript Fig.3)
#| 
#| # 3. FIGURE 2 — BETA DIVERSITY
#| # WEIGHTED + UNWEIGHTED UNIFRAC PCoA
#| # ============================================================
#| 
#| set.seed(123)
#| 
#| ord_w <- ordinate(
#|   prok_data,
#|   method = "PCoA",
#|   distance = "wunifrac"
#| )
#| 
#| ord_u <- ordinate(
#|   prok_data,
#|   method = "PCoA",
#|   distance = "unifrac"
#| )
#| 
#| # Function to extract PCoA coordinates and variance
#| get_pcoa <- function(ord, meta) {
#|   
#|   dat <- as.data.frame(
#|     ord$vectors[,1:2]
#|   ) %>%
#|     rownames_to_column("Sample") %>%
#|     left_join( meta %>% select( Sample, treatment, time, electrode ), by = "Sample" )
#|   
#|   if (!is.null(ord$values$Relative_eig)) {
#|     
#|     var_exp <- ord$values$Relative_eig[1:2] * 100
#|     
#|   } else {
#|     
#|     eig <- ord$values$Eigenvalues
#|     
#|     var_exp <- eig[1:2] /
#|       sum(eig[eig > 0]) * 100
#|   }
#|   
#|   list(
#|     data = dat,
#|     var = var_exp
#|   )
#| }
#| 
#| pcoa_w <- get_pcoa(
#|   ord_w,
#|   meta
#| )
#| 
#| pcoa_u <- get_pcoa(
#|   ord_u,
#|   meta
#| )
#| 
#| # H/O pairs
#| connectors <- meta %>%
#|   filter( electrode %in% c("H","O") ) %>%
#|   select( Sample, treatment, time, electrode )
#| 
#| make_pcoa_plot <- function(obj, title) {
#|   
#|   dat <- obj$data
#|   
#|   seg <- connectors %>%
#|     left_join( dat %>% select( Sample, Axis.1, Axis.2 ), by = "Sample" ) %>%
#|     pivot_wider( names_from = electrode, values_from = c( Sample, Axis.1, Axis.2 ) ) %>%
#|     filter( !is.na(Axis.1_H), !is.na(Axis.1_O) )
#|   
#|   ggplot(
#|     dat,
#|     aes( x = Axis.1, y = Axis.2, colour = treatment, shape = electrode )
#|   ) +
#|     
#|     geom_segment( data = seg, aes( x = Axis.1_H, y = Axis.2_H, xend = Axis.1_O, yend = Axis.2_O ), inherit.aes = FALSE, linetype = 2, linewidth = .5 ) +
#|     
#|     geom_point( size = 4 ) +
#|     
#|     geom_text( aes( label = Sample ), vjust = -1, size = 3, show.legend = FALSE ) +
#|     
#|     labs( title = title, x = paste0( "PCoA1 (", round(obj$var[1],1), "%)" ), y = paste0( "PCoA2 (", round(obj$var[2],1), "%)" ) ) +
#|     
#|     theme_classic( base_size = 11 )
#| }
#| 
#| p_w <- make_pcoa_plot(
#|   pcoa_w,
#|   "Weighted UniFrac"
#| )
#| 
#| p_u <- make_pcoa_plot(
#|   pcoa_u,
#|   "Unweighted UniFrac"
#| )
#| 
#| p_beta <- p_w + p_u +
#|   plot_layout(
#|     guides = "collect"
#|   ) &
#|   theme( legend.position = "bottom" )
#| 
#| p_beta
#| 
#| save_fig(
#|   p_beta,
#|   "Figure2_Beta_diversity_UniFrac",
#|   11,
#|   5.5
#| )
#| 
#| # ============================================================
#| # 3B. H–O UNIFRAC DISTANCES — RESULTS TABLE
#| # ============================================================
#| 
#| dist_w <- phyloseq::distance(
#|   prok_data,
#|   method = "wunifrac"
#| )
#| 
#| dist_u <- phyloseq::distance(
#|   prok_data,
#|   method = "unifrac"
#| )
#| 
#| dist_w_mat <- as.matrix(dist_w)
#| dist_u_mat <- as.matrix(dist_u)
#| 
#| HO_pairs <- meta %>%
#|   filter( electrode %in% c("H","O") ) %>%
#|   select( Sample, treatment, time, electrode ) %>%
#|   pivot_wider( names_from = electrode, values_from = Sample ) %>%
#|   filter( !is.na(H), !is.na(O) )
#| 
#| HO_distances <- HO_pairs %>%
#|   rowwise() %>%
#|   mutate( Weighted_UniFrac = dist_w_mat[H,O], Unweighted_UniFrac = dist_u_mat[H,O] ) %>%
#|   ungroup()
#| 
#| cat("\n================ H-O UNIFRAC DISTANCES ================\n")
#| 
#| print( HO_distances, n = Inf )
#| 
#| write.csv( HO_distances, "Table_HO_UniFrac_distances.csv", row.names = FALSE )
#| 
#| # ============================================================
#| # 3C. EXPLORATORY STATISTICS — UNREPLICATED DESIGN
#| # ============================================================
#| #
#| # IMPORTANT INTERPRETATIVE NOTE
#| # The experimental units were not independently replicated.
#| # Therefore, inferential tests must not be interpreted as
#| # confirmatory tests of treatment effects. The analyses below
#| # are intended to quantify effect sizes, temporal trajectories,
#| # multivariate separation and dispersion descriptively/
#| # exploratorily.
#| #
#| # In particular:
#| # - no factorial Treatment x Time x Electrode PERMANOVA is used;
#| # - H and O are compartments of the same experimental system,
#| #   not independent biological replicates;
#| # - V1OT1 (PCR-negative / no usable post-filtering sequences)
#| #   must remain absent rather than being coded as zero diversity.
#| # ============================================================
#| 
#| dir.create("Exploratory_statistics", showWarnings = FALSE)
#| 
#| # ------------------------------------------------------------
#| # A. Experimental-design audit
#| # ------------------------------------------------------------
#| 
#| design_audit <- meta %>%

# ===== ORIGINAL SOURCE: 05_Fig4_tassonomia_originale_finale_v4.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 15820-16967; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Order composition and Bray-Curtis dendrogram (original Fig.3 -> manuscript Fig.4)
#| 
#| # FIGURE 3 - VERSION 4
#| #
#| # ORDER-LEVEL PROKARYOTIC COMMUNITY COMPOSITION
#| # AND BRAY-CURTIS UPGMA CLUSTERING
#| #
#| # FINAL CLEAN VERSION FOR COMPARISON
#| #
#| # ANALYTICAL STRATEGY:
#| #   - Bray-Curtis calculated using ALL Orders
#| #   - UPGMA hierarchical clustering (average linkage)
#| #   - Orders reaching >= 5% in at least one sample shown individually
#| #   - Remaining Orders grouped as "Other" ONLY for visualization
#| #
#| # GRAPHICAL STRATEGY:
#| #   - exact stacked-bar coordinates calculated manually
#| #   - numerical codes assigned to dominant Orders
#| #   - numerical codes shown only for segments >= 3%
#| #   - numerical codes centred geometrically within segments
#| #   - all numerical codes shown in white bold text
#| #   - slightly darker colours for light taxa
#| #   - "Other" shown in medium grey and not numbered
#| #   - complete 0-100% bar outlined to make terminal segments visible
#| #   - legend positioned below the figure
#| #
#| # INPUT:
#| #   order_all
#| #
#| # OUTPUT:
#| #   Fig3_Order_BrayCurtis_UPGMA_numbered_V4.pdf
#| #   Fig3_Order_BrayCurtis_UPGMA_numbered_V4.svg
#| #   Fig3_Order_BrayCurtis_UPGMA_numbered_V4.png
#| ###############################################################################
#| 
#| 
#| # =============================================================================
#| # 0. PACKAGES
#| # =============================================================================
#| 
#| library(tidyverse)
#| library(vegan)
#| library(ggdendro)
#| library(patchwork)
#| library(grid)
#| 
#| 
#| # =============================================================================
#| # 1. PREPARE ORDER-LEVEL COMMUNITY MATRIX
#| # =============================================================================
#| 
#| order_mat_fig3_v4 <- as.data.frame(order_all)
#| 
#| # Retain numeric abundance columns only
#| order_mat_fig3_v4 <- order_mat_fig3_v4[
#|   ,
#|   sapply(order_mat_fig3_v4, is.numeric),
#|   drop = FALSE
#| ]
#| 
#| # Remove Orders absent from all samples
#| order_mat_fig3_v4 <- order_mat_fig3_v4[
#|   ,
#|   colSums(order_mat_fig3_v4, na.rm = TRUE) > 0,
#|   drop = FALSE
#| ]
#| 
#| order_mat_fig3_v4 <- as.matrix(order_mat_fig3_v4)
#| 
#| storage.mode(order_mat_fig3_v4) <- "numeric"
#| 
#| 
#| # =============================================================================
#| # 2. CHECK SAMPLE TOTALS
#| # =============================================================================
#| 
#| row_totals_fig3_v4 <- rowSums(
#|   order_mat_fig3_v4,
#|   na.rm = TRUE
#| )
#| 
#| if (any(row_totals_fig3_v4 == 0)) {
#|   
#|   stop(
#|     "At least one sample has total Order abundance equal to zero."
#|   )
#| }
#| 
#| 
#| # =============================================================================
#| # 3. CONVERT TO RELATIVE ABUNDANCE
#| # =============================================================================
#| 
#| order_rel_fig3_v4 <- sweep(
#|   order_mat_fig3_v4,
#|   1,
#|   row_totals_fig3_v4,
#|   "/"
#| )
#| 
#| 
#| cat("\n============================================================\n")
#| cat("FIGURE 3 V4 - RELATIVE-ABUNDANCE ROW SUMS\n")
#| cat("============================================================\n")
#| 
#| print(
#|   rowSums(order_rel_fig3_v4)
#| )
#| 
#| 
#| # =============================================================================
#| # 4. BRAY-CURTIS DISSIMILARITY
#| #
#| # IMPORTANT:
#| # ALL retained Orders are used here.
#| # No abundance threshold is applied before Bray-Curtis calculation.
#| # =============================================================================
#| 
#| bray_fig3_v4 <- vegdist(
#|   order_rel_fig3_v4,
#|   method = "bray"
#| )
#| 
#| 
#| cat("\n============================================================\n")
#| cat("FIGURE 3 V4 - BRAY-CURTIS MATRIX\n")
#| cat("============================================================\n")
#| 
#| print(
#|   round(
#|     as.matrix(bray_fig3_v4),
#|     3
#|   )
#| )
#| 
#| 
#| # =============================================================================
#| # 5. UPGMA HIERARCHICAL CLUSTERING
#| # =============================================================================
#| 
#| clust_fig3_v4 <- hclust(
#|   bray_fig3_v4,
#|   method = "average"
#| )
#| 
#| 
#| dend_fig3_v4 <- as.dendrogram(
#|   clust_fig3_v4
#| )
#| 
#| 
#| ddata_fig3_v4 <- dendro_data(
#|   dend_fig3_v4,
#|   type = "rectangle"
#| )
#| 
#| 
#| # =============================================================================
#| # 6. SAMPLE ORDER FROM DENDROGRAM
#| # =============================================================================
#| 
#| sample_order_fig3_v4 <- ddata_fig3_v4$labels %>%
#|   
#|   arrange(x) %>%
#|   
#|   pull(label)
#| 
#| 
#| cat("\n============================================================\n")
#| cat("FIGURE 3 V4 - DENDROGRAM SAMPLE ORDER\n")
#| cat("============================================================\n")
#| 
#| print(
#|   sample_order_fig3_v4
#| )
#| 
#| 
#| # =============================================================================
#| # 7. SELECT DOMINANT ORDERS FOR VISUALIZATION
#| #
#| # Criterion:
#| # >= 5% relative abundance in at least one sample
#| #
#| # IMPORTANT:
#| # This threshold affects ONLY graphical representation.
#| # =============================================================================
#| 
#| order_max_fig3_v4 <- apply(
#|   order_rel_fig3_v4,
#|   2,
#|   max,
#|   na.rm = TRUE
#| )
#| 
#| 
#| dominant_orders_fig3_v4 <- names(
#|   order_max_fig3_v4[
#|     order_max_fig3_v4 >= 0.05
#|   ]
#| )
#| 
#| 
#| # Alphabetical order gives stable and reproducible numbering
#| dominant_orders_fig3_v4 <- sort(
#|   dominant_orders_fig3_v4
#| )
#| 
#| 
#| cat("\n============================================================\n")
#| cat("FIGURE 3 V4 - DOMINANT ORDERS >= 5%\n")
#| cat("============================================================\n")
#| 
#| print(
#|   dominant_orders_fig3_v4
#| )
#| 
#| cat(
#|   "\nNumber of dominant Orders:",
#|   length(dominant_orders_fig3_v4),
#|   "\n"
#| )
#| 
#| 
#| # =============================================================================
#| # 8. ASSIGN NUMERICAL TAXONOMIC CODES
#| #
#| # IMPORTANT:
#| # Codes are generated HERE from the same sorted vector used later
#| # for colours, factor levels and legend.
#| #
#| # "Other" is deliberately not numbered.
#| # =============================================================================
#| 
#| order_codes_fig3_v4 <- tibble(
#|   
#|   Order_plot = dominant_orders_fig3_v4,
#|   
#|   Taxon_number = seq_along(
#|     dominant_orders_fig3_v4
#|   )
#|   
#| ) %>%
#|   
#|   mutate(
#|     
#|     Legend_label = paste0(
#|       Taxon_number,
#|       "  ",
#|       Order_plot
#|     )
#|   )
#| 
#| 
#| cat("\n============================================================\n")
#| cat("FIGURE 3 V4 - TAXONOMIC CODES\n")
#| cat("============================================================\n")
#| 
#| print(
#|   order_codes_fig3_v4,
#|   n = Inf
#| )
#| 
#| 
#| # =============================================================================
#| # 9. PREPARE STACKED-BAR DATA
#| # =============================================================================
#| 
#| bar_data_fig3_v4 <- as.data.frame(
#|   order_rel_fig3_v4
#| ) %>%
#|   
#|   rownames_to_column(
#|     "Sample"
#|   ) %>%
#|   
#|   pivot_longer(
#|     
#|     cols = -Sample,
#|     
#|     names_to = "Order",
#|     
#|     values_to = "Relative_abundance"
#|   ) %>%
#|   
#|   mutate(
#|     
#|     Order_plot = if_else(
#|       
#|       Order %in% dominant_orders_fig3_v4,
#|       
#|       Order,
#|       
#|       "Other"
#|     )
#|   ) %>%
#|   
#|   group_by(
#|     Sample,
#|     Order_plot
#|   ) %>%
#|   
#|   summarise(
#|     
#|     Relative_abundance = sum(
#|       Relative_abundance,
#|       na.rm = TRUE
#|     ),
#|     
#|     .groups = "drop"
#|   ) %>%
#|   
#|   mutate(
#|     
#|     Relative_abundance =
#|       Relative_abundance * 100
#|   )
#| 
#| 
#| # =============================================================================
#| # 10. ADD NUMERICAL TAXONOMIC CODES
#| # =============================================================================
#| 
#| bar_data_fig3_v4 <- bar_data_fig3_v4 %>%
#|   
#|   left_join(
#|     
#|     order_codes_fig3_v4 %>%
#|       
#|       select(
#|         Order_plot,
#|         Taxon_number
#|       ),
#|     
#|     by = "Order_plot"
#|   )
#| 
#| 
#| # =============================================================================
#| # 11. CHECK THAT EVERY SAMPLE SUMS TO 100%
#| # =============================================================================
#| 
#| fig3_totals_v4 <- bar_data_fig3_v4 %>%
#|   
#|   group_by(Sample) %>%
#|   
#|   summarise(
#|     
#|     Total = sum(
#|       Relative_abundance,
#|       na.rm = TRUE
#|     ),
#|     
#|     .groups = "drop"
#|   )
#| 
#| 
#| cat("\n============================================================\n")
#| cat("FIGURE 3 V4 - BAR TOTALS\n")
#| cat("============================================================\n")
#| 
#| fig3_totals_v4 %>%
#|   
#|   mutate(
#|     
#|     Total_print =
#|       sprintf(
#|         "%.10f",
#|         Total
#|       )
#|   ) %>%
#|   
#|   print(
#|     n = Inf
#|   )
#| 
#| 
#| if (
#|   any(
#|     abs(fig3_totals_v4$Total - 100) > 1e-6
#|   )
#| ) {
#|   
#|   stop(
#|     "ERROR: one or more samples do not sum to 100%."
#|   )
#| }
#| 
#| 
#| # =============================================================================
#| # 12. DEFINE TAXONOMIC STACKING ORDER
#| #
#| # Dominant Orders first, Other last.
#| # =============================================================================
#| 
#| order_levels_fig3_v4 <- c(
#|   
#|   dominant_orders_fig3_v4,
#|   
#|   "Other"
#| )
#| 
#| 
#| bar_data_fig3_v4 <- bar_data_fig3_v4 %>%
#|   
#|   mutate(
#|     
#|     Order_plot = factor(
#|       Order_plot,
#|       levels = order_levels_fig3_v4
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 13. CALCULATE EXACT STACKED-SEGMENT COORDINATES
#| #
#| # Coordinates are calculated on the COMPLETE bar before label filtering.
#| #
#| # xmin = beginning of segment
#| # xmax = end of segment
#| # xmid = exact geometric centre of segment
#| # =============================================================================
#| 
#| bar_data_fig3_v4 <- bar_data_fig3_v4 %>%
#|   
#|   arrange(
#|     Sample,
#|     Order_plot
#|   ) %>%
#|   
#|   group_by(Sample) %>%
#|   
#|   mutate(
#|     
#|     xmin =
#|       cumsum(Relative_abundance) -
#|       Relative_abundance,
#|     
#|     xmax =
#|       cumsum(Relative_abundance),
#|     
#|     xmid =
#|       (
#|         xmin +
#|           xmax
#|       ) / 2
#|   ) %>%
#|   
#|   ungroup()
#| 
#| 
#| # =============================================================================
#| # 14. VERIFY FINAL XMAX
#| #
#| # Every sample MUST terminate at 100.
#| # =============================================================================
#| 
#| fig3_coordinates_check_v4 <- bar_data_fig3_v4 %>%
#|   
#|   group_by(Sample) %>%
#|   
#|   summarise(
#|     
#|     Total =
#|       sum(Relative_abundance),
#|     
#|     Final_xmax =
#|       max(xmax),
#|     
#|     .groups = "drop"
#|   )
#| 
#| 
#| cat("\n============================================================\n")
#| cat("FIGURE 3 V4 - FINAL COORDINATE CHECK\n")
#| cat("============================================================\n")
#| 
#| print(
#|   fig3_coordinates_check_v4,
#|   n = Inf
#| )
#| 
#| 
#| if (
#|   any(
#|     abs(fig3_coordinates_check_v4$Final_xmax - 100) > 1e-6
#|   )
#| ) {
#|   
#|   stop(
#|     "ERROR: one or more graphical bars do not terminate at x = 100."
#|   )
#| }
#| 
#| 
#| # =============================================================================
#| # 15. EXACT VERTICAL POSITIONS FROM DENDROGRAM
#| #
#| # Same y coordinates are used for dendrogram leaves and composition bars.
#| # =============================================================================
#| 
#| sample_positions_fig3_v4 <- ddata_fig3_v4$labels %>%
#|   
#|   transmute(
#|     
#|     Sample = as.character(label),
#|     
#|     ypos = x
#|   )
#| 
#| 
#| bar_rect_fig3_v4 <- bar_data_fig3_v4 %>%
#|   
#|   mutate(
#|     Sample = as.character(Sample)
#|   ) %>%
#|   
#|   left_join(
#|     
#|     sample_positions_fig3_v4,
#|     
#|     by = "Sample"
#|   )
#| 
#| 
#| # Safety check
#| if (
#|   any(
#|     is.na(bar_rect_fig3_v4$ypos)
#|   )
#| ) {
#|   
#|   stop(
#|     "ERROR: sample names in bar data and dendrogram do not match."
#|   )
#| }
#| 
#| 
#| # =============================================================================
#| # 16. LABEL DATA
#| #
#| # Numerical code displayed only for segments >= 3%.
#| #
#| # This is exclusively a graphical readability threshold.
#| # It has no analytical meaning.
#| #
#| # Other is never numbered.
#| # =============================================================================
#| 
#| label_data_fig3_v4 <- bar_rect_fig3_v4 %>%
#|   
#|   filter(
#|     
#|     as.character(Order_plot) != "Other",
#|     
#|     Relative_abundance >= 3
#|   )
#| 
#| 
#| # =============================================================================
#| # 17. HIGH-CONTRAST PALETTE
#| #
#| # V4 modifications:
#| #
#| # - Taxon 7 is darker than in V3
#| # - Taxon 16 is darker than in V3
#| # - Other = grey70
#| #
#| # This improves visibility of white numerical labels.
#| # =============================================================================
#| 
#| high_contrast_colors_v4 <- c(
#|   
#|   "#E41A1C",  # 1
#|   "#377EB8",  # 2
#|   "#4DAF4A",  # 3
#|   "#984EA3",  # 4
#|   "#FF7F00",  # 5
#|   
#|   "#A65628",  # 6
#|   
#|   "#D94F9D",  # 7 - darker pink
#|   
#|   "#1B9E77",  # 8
#|   "#D95F02",  # 9
#|   "#7570B3",  # 10
#|   
#|   "#66A61E",  # 11
#|   "#E7298A",  # 12
#|   "#A6761D",  # 13
#|   "#1F78B4",  # 14
#|   "#33A02C",  # 15
#|   
#|   "#E76F73",  # 16 - darker salmon
#|   
#|   "#6A3D9A",  # 17
#|   "#B15928",  # 18
#|   "#17BECF",  # 19
#|   "#BCBD22"   # 20
#| )
#| 
#| 
#| if (
#|   length(dominant_orders_fig3_v4) >
#|   length(high_contrast_colors_v4)
#| ) {
#|   
#|   stop(
#|     paste(
#|       "More dominant Orders than available colours.",
#|       "Extend high_contrast_colors_v4."
#|     )
#|   )
#| }
#| 
#| 
#| tax_colors_fig3_v4 <- c(
#|   
#|   setNames(
#|     
#|     high_contrast_colors_v4[
#|       seq_along(
#|         dominant_orders_fig3_v4
#|       )
#|     ],
#|     
#|     dominant_orders_fig3_v4
#|   ),
#|   
#|   Other = "grey70"
#| )
#| 
#| 
#| # =============================================================================
#| # 18. NUMBERED LEGEND LABELS
#| #
#| # Generated directly from order_codes_fig3_v4.
#| # =============================================================================
#| 
#| legend_labels_fig3_v4 <- c(
#|   
#|   setNames(
#|     
#|     order_codes_fig3_v4$Legend_label,
#|     
#|     order_codes_fig3_v4$Order_plot
#|   ),
#|   
#|   Other = "Other"
#| )
#| 
#| 
#| # =============================================================================
#| # 19. FINAL TAXON NUMBER / COLOUR / LEGEND CHECK
#| #
#| # Useful for tracking and for avoiding inconsistencies between R objects.
#| # =============================================================================
#| 
#| fig3_taxon_key_v4 <- order_codes_fig3_v4 %>%
#|   
#|   mutate(
#|     
#|     Colour =
#|       tax_colors_fig3_v4[
#|         Order_plot
#|       ]
#|   )
#| 
#| 
#| cat("\n============================================================\n")
#| cat("FIGURE 3 V4 - FINAL TAXON KEY\n")
#| cat("============================================================\n")
#| 
#| print(
#|   fig3_taxon_key_v4,
#|   n = Inf
#| )
#| 
#| 
#| # =============================================================================
#| # 20. LEFT PANEL
#| # HORIZONTAL BRAY-CURTIS UPGMA DENDROGRAM
#| # =============================================================================
#| 
#| max_dist_fig3_v4 <- max(
#|   ddata_fig3_v4$segments$y,
#|   na.rm = TRUE
#| )
#| 
#| 
#| p_dend_fig3_v4 <- ggplot() +
#|   
#|   # ---------------------------------------------------------------------------
#| # Dendrogram branches
#| # ---------------------------------------------------------------------------
#| 
#| geom_segment(
#|   
#|   data = ddata_fig3_v4$segments,
#|   
#|   aes(
#|     x = y,
#|     y = x,
#|     xend = yend,
#|     yend = xend
#|   ),
#|   
#|   linewidth = 0.55,
#|   
#|   lineend = "square"
#| ) +
#|   
#|   # ---------------------------------------------------------------------------
#| # Sample labels
#| # ---------------------------------------------------------------------------
#| 
#| geom_text(
#|   
#|   data = ddata_fig3_v4$labels,
#|   
#|   aes(
#|     x = -0.055,
#|     y = x,
#|     label = label
#|   ),
#|   
#|   hjust = 0,
#|   
#|   size = 3.2,
#|   
#|   fontface = "plain"
#| ) +
#|   
#|   # ---------------------------------------------------------------------------
#| # Bray-Curtis axis
#| # ---------------------------------------------------------------------------
#| 
#| scale_x_reverse(
#|   
#|   limits = c(
#|     max_dist_fig3_v4,
#|     -0.20
#|   ),
#|   
#|   expand = c(
#|     0,
#|     0
#|   )
#| ) +
#|   
#|   # ---------------------------------------------------------------------------
#| # Exact vertical coordinates
#| # ---------------------------------------------------------------------------
#| 
#| scale_y_continuous(
#|   
#|   limits = c(
#|     0.5,
#|     length(sample_order_fig3_v4) + 0.5
#|   ),
#|   
#|   breaks = NULL,
#|   
#|   expand = c(
#|     0,
#|     0
#|   )
#| ) +
#|   
#|   coord_cartesian(
#|     clip = "off"
#|   ) +
#|   
#|   labs(
#|     
#|     x = "Bray-Curtis dissimilarity",
#|     
#|     y = NULL
#|   ) +
#|   
#|   theme_classic(
#|     base_size = 10
#|   ) +
#|   
#|   theme(
#|     
#|     axis.line.y =
#|       element_blank(),
#|     
#|     axis.ticks.y =
#|       element_blank(),
#|     
#|     axis.text.x =
#|       element_text(
#|         colour = "black",
#|         size = 8
#|       ),
#|     
#|     axis.title.x =
#|       element_text(
#|         size = 9
#|       ),
#|     
#|     plot.margin =
#|       margin(
#|         t = 5,
#|         r = 35,
#|         b = 5,
#|         l = 5
#|       )
#|   )
#| 
#| 
#| # =============================================================================
#| # =============================================================================
#| # 21. RIGHT PANEL
#| # EXACT HORIZONTAL STACKED BARPLOT
#| #
#| # ROBUST STRATEGY:
#| #   1. draw the complete 0-100% bar in grey ("Other")
#| #   2. overlay all individually represented Orders
#| #   3. leave the residual fraction visible as grey
#| #
#| # This prevents any terminal white gap caused by rendering.
#| # =============================================================================
#| 
#| 
#| # Data containing ONLY individually represented Orders
#| dominant_rect_fig3_v4 <- bar_rect_fig3_v4 %>%
#|   
#|   filter(
#|     as.character(Order_plot) != "Other",
#|     Relative_abundance > 0
#|   )
#| 
#| 
#| p_bar_fig3_v4 <- ggplot() +
#|   
#|   # ---------------------------------------------------------------------------
#| # BACKGROUND:
#| # complete 0-100% bar = Other
#| #
#| # Every sample is first drawn as a continuous grey bar.
#| # ---------------------------------------------------------------------------
#| 
#| geom_rect(
#|   
#|   data = sample_positions_fig3_v4,
#|   
#|   aes(
#|     xmin = 0,
#|     xmax = 100,
#|     ymin = ypos - 0.38,
#|     ymax = ypos + 0.38
#|   ),
#|   
#|   inherit.aes = FALSE,
#|   
#|   fill = "grey70",
#|   
#|   colour = NA
#| ) +
#|   
#|   # ---------------------------------------------------------------------------
#| # DOMINANT ORDERS
#| #
#| # Drawn over the grey background using their exact coordinates.
#| # ---------------------------------------------------------------------------
#| 
#| geom_rect(
#|   
#|   data = dominant_rect_fig3_v4,
#|   
#|   aes(
#|     xmin = xmin,
#|     xmax = xmax,
#|     ymin = ypos - 0.38,
#|     ymax = ypos + 0.38,
#|     fill = Order_plot
#|   ),
#|   
#|   colour = "white",
#|   
#|   linewidth = 0.15
#| ) +
#|   
#|   # ---------------------------------------------------------------------------
#| # COMPLETE BAR OUTLINE
#| # ---------------------------------------------------------------------------
#| 
#| geom_rect(
#|   
#|   data = sample_positions_fig3_v4,
#|   
#|   aes(
#|     xmin = 0,
#|     xmax = 100,
#|     ymin = ypos - 0.38,
#|     ymax = ypos + 0.38
#|   ),
#|   
#|   inherit.aes = FALSE,
#|   
#|   fill = NA,
#|   
#|   colour = "grey35",
#|   
#|   linewidth = 0.25
#| ) +
#|   
#|   # ---------------------------------------------------------------------------
#| # NUMERICAL TAXONOMIC CODES
#| # ---------------------------------------------------------------------------
#| 
#| geom_text(
#|   
#|   data = label_data_fig3_v4,
#|   
#|   aes(
#|     x = xmid,
#|     y = ypos,
#|     label = Taxon_number
#|   ),
#|   
#|   inherit.aes = FALSE,
#|   
#|   colour = "white",
#|   
#|   fontface = "bold",
#|   
#|   size = 2.8,
#|   
#|   hjust = 0.5,
#|   
#|   vjust = 0.5,
#|   
#|   show.legend = FALSE
#| ) +
#|   
#|   # ---------------------------------------------------------------------------
#| # TAXONOMIC COLOURS
#| # ---------------------------------------------------------------------------
#| 
#| scale_fill_manual(
#|   
#|   values = tax_colors_fig3_v4,
#|   
#|   breaks = order_levels_fig3_v4,
#|   
#|   labels = legend_labels_fig3_v4,
#|   
#|   drop = FALSE
#| ) +
#|   
#|   # ---------------------------------------------------------------------------
#| # X AXIS
#| # ---------------------------------------------------------------------------
#| 
#| scale_x_continuous(
#|   
#|   limits = c(0, 100),
#|   
#|   breaks = c(
#|     0,
#|     25,
#|     50,
#|     75,
#|     100
#|   ),
#|   
#|   expand = c(0, 0)
#| ) +
#|   
#|   # ---------------------------------------------------------------------------
#| # Y AXIS
#| # ---------------------------------------------------------------------------
#| 
#| scale_y_continuous(
#|   
#|   limits = c(
#|     0.5,
#|     length(sample_order_fig3_v4) + 0.5
#|   ),
#|   
#|   breaks = NULL,
#|   
#|   expand = c(0, 0)
#| ) +
#|   
#|   labs(
#|     
#|     x = "Relative abundance (%)",
#|     
#|     y = NULL,
#|     
#|     fill = "Order"
#|   ) +
#|   
#|   theme_classic(
#|     base_size = 10
#|   ) +
#|   
#|   theme(
#|     
#|     axis.line.y = element_blank(),
#|     
#|     axis.ticks.y = element_blank(),
#|     
#|     axis.text.y = element_blank(),
#|     
#|     axis.text.x = element_text(
#|       colour = "black",
#|       size = 8
#|     ),
#|     
#|     axis.title.x = element_text(
#|       size = 9
#|     ),
#|     
#|     plot.margin = margin(
#|       t = 5,
#|       r = 5,
#|       b = 5,
#|       l = 0
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 22. COMBINE DENDROGRAM + TAXONOMIC COMPOSITION
#| # =============================================================================
#| 
#| fig3_v4 <- p_dend_fig3_v4 +
#|   
#|   p_bar_fig3_v4 +
#|   
#|   plot_layout(
#|     
#|     widths = c(
#|       1.30,
#|       3.70
#|     ),
#|     
#|     guides = "collect"
#|   ) &
#|   
#|   theme(
#|     
#|     legend.position =
#|       "bottom",
#|     
#|     legend.title =
#|       element_text(
#|         face = "bold",
#|         size = 9
#|       ),
#|     
#|     legend.text =
#|       element_text(
#|         size = 8
#|       ),
#|     
#|     legend.key.width =
#|       unit(
#|         0.40,
#|         "cm"
#|       ),
#|     
#|     legend.key.height =
#|       unit(
#|         0.35,
#|         "cm"
#|       )
#|   ) &
#|   
#|   guides(
#|     
#|     fill = guide_legend(
#|       
#|       nrow = 3,
#|       
#|       byrow = TRUE,
#|       
#|       title.position = "top"
#|     )
#|   )
#| 
#| 
#| # =============================================================================
#| # 23. DISPLAY FINAL FIGURE
#| # =============================================================================
#| 
#| print(
#|   fig3_v4
#| )
#| 
#| 
#| # =============================================================================
#| # 24. EXPORT
#| # =============================================================================
#| 
#| ggsave(
#|   
#|   filename =
#|     "Fig3_Order_BrayCurtis_UPGMA_numbered_V4.pdf",
#|   
#|   plot = fig3_v4,
#|   
#|   width = 12,
#|   
#|   height = 7.4,
#|   
#|   units = "in"
#| )
#| 
#| 
#| ggsave(
#|   
#|   filename =
#|     "Fig3_Order_BrayCurtis_UPGMA_numbered_V4.svg",
#|   
#|   plot = fig3_v4,
#|   
#|   width = 12,
#|   
#|   height = 7.4,
#|   
#|   units = "in"
#| )
#| 
#| 
#| ggsave(
#|   
#|   filename =
#|     "Fig3_Order_BrayCurtis_UPGMA_numbered_V4.png",
#|   
#|   plot = fig3_v4,
#|   
#|   width = 12,
#|   
#|   height = 7.4,
#|   
#|   units = "in",
#|   
#|   dpi = 600
#| )
#| 
#| 
#| ###############################################################################
#| # END FIGURE 3 - VERSION 4

# ===== ORIGINAL SOURCE: 06_Fig5_heatmap_ordini.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 6969-7053; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Order z-score heatmap (original Fig.3B -> manuscript Fig.5)
#| 
#| # 6. FIGURE 3B — ORDER × SAMPLE Z-SCORE HEATMAP
#| # ============================================================
#| 
#| # More permissive threshold:
#| # retain condition-specific Orders reaching >=0.5%
#| orders_keep_heat <- order_abund %>%
#|   filter( Order != "Unclassified" ) %>%
#|   group_by(Order) %>%
#|   summarise( Max = max(Abundance), .groups = "drop" ) %>%
#|   filter( Max >= 0.5 ) %>%
#|   pull(Order)
#| 
#| order_heat <- order_abund %>%
#|   filter( Order != "Unclassified", Order %in% orders_keep_heat ) %>%
#|   complete(
#|     Sample = sample_names(prok_data),
#|     Order = orders_keep_heat,
#|     fill = list(
#|       Abundance = 0
#|     )
#|   ) %>%
#|   group_by(Order) %>%
#|   mutate( Z = ( Abundance - mean(Abundance) ) / sd(Abundance) ) %>%
#|   ungroup() %>%
#|   filter( is.finite(Z) ) %>%
#|   left_join( meta %>% select( Sample, treatment, time, electrode ), by = "Sample" )
#| 
#| order_heat$Sample <- factor(
#|   order_heat$Sample,
#|   levels = sample_levels
#| )
#| 
#| # Cluster Orders based on their z-score profiles
#| z_matrix <- order_heat %>%
#|   select( Sample, Order, Z ) %>%
#|   pivot_wider( names_from = Sample, values_from = Z, values_fill = 0 ) %>%
#|   column_to_rownames("Order") %>%
#|   as.matrix()
#| 
#| hc <- hclust(
#|   dist(z_matrix),
#|   method = "average"
#| )
#| 
#| order_levels_heat <- hc$labels[
#|   hc$order
#| ]
#| 
#| order_heat$Order <- factor(
#|   order_heat$Order,
#|   levels = rev(
#|     order_levels_heat
#|   )
#| )
#| 
#| p_order_heat <- ggplot(
#|   order_heat,
#|   aes( x = Sample, y = Order, fill = Z )
#| ) +
#|   
#|   geom_tile( colour = "white", linewidth = .25 ) +
#|   
#|   scale_fill_gradient2( low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-3,3), oob = squish, name = "Z-score" ) +
#|   
#|   labs( x = NULL, y = NULL ) +
#|   
#|   theme_classic( base_size = 10 ) +
#|   
#|   theme(
#|     axis.text.x = element_text(
#|       angle = 45,
#|       hjust = 1,
#|       face = "bold",
#|       colour = "black"
#|     ),
#|     axis.text.y = element_text(
#|       size = 7.5,
#|       colour = "black"
#|     ),
#|     axis.ticks = element_blank()
#|   )
#| 
#| p_order_heat
#| 
#| # ============================================================

# ===== ORIGINAL SOURCE: 07_Fig6_nMDS_envfit_dbRDA.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 7188-7410; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Bray-Curtis NMDS envfit univariate dbRDA (original Fig.4 -> manuscript Fig.6)
#| 
#| # 10. FIGURE 4A — BRAY-CURTIS nMDS + ENVFIT
#| # ============================================================
#| 
#| set.seed(123)
#| 
#| bac_nmds <- metaMDS(
#|   bac_ra.v,
#|   distance = "bray",
#|   trymax = 100,
#|   autotransform = FALSE
#| )
#| 
#| nmds_scores <- as.data.frame(
#|   scores(
#|     bac_nmds,
#|     display = "sites"
#|   )
#| ) %>%
#|   rownames_to_column(
#|     "Sample"
#|   ) %>%
#|   left_join( meta %>% select( Sample, treatment, time, electrode ), by = "Sample" )
#| 
#| # ENVFIT
#| env_fit <- envfit(
#|   bac_nmds,
#|   prok_sample[
#|     rownames(bac_ra.v),
#|     env_vars
#|   ],
#|   permutations = 999
#| )
#| 
#| vec <- as.data.frame(
#|   scores(
#|     env_fit,
#|     display = "vectors"
#|   )
#| ) %>%
#|   rownames_to_column(
#|     "Variable"
#|   )
#| 
#| vec$R2 <- env_fit$vectors$r
#| vec$p <- env_fit$vectors$pvals
#| vec$FDR <- p.adjust(
#|   vec$p,
#|   method = "BH"
#| )
#| 
#| cat("\n================ ENVFIT ================\n")
#| 
#| vec %>%
#|   as_tibble() %>%
#|   print(n = Inf)
#| 
#| write.csv( vec, "Table_envfit.csv", row.names = FALSE )
#| 
#| # Show nominal p < 0.05 vectors
#| # FDR significance is reported separately
#| vec_plot <- vec %>%
#|   filter( p < .05 )
#| 
#| if (nrow(vec_plot) > 0) {
#|   
#|   mult <- vegan:::ordiArrowMul(
#|     as.matrix(
#|       vec_plot[
#|         ,
#|         c(
#|           "NMDS1",
#|           "NMDS2"
#|         )
#|       ]
#|     )
#|   )
#|   
#|   vec_plot <- vec_plot %>%
#|     mutate( NMDS1 = NMDS1 * mult, NMDS2 = NMDS2 * mult )
#| }
#| 
#| # H-O connectors
#| seg_nmds <- meta %>%
#|   
#|   filter( electrode %in% c( "H", "O" ) ) %>%
#|   
#|   select( Sample, treatment, time, electrode ) %>%
#|   
#|   left_join( nmds_scores %>% select( Sample, NMDS1, NMDS2 ), by = "Sample" ) %>%
#|   
#|   pivot_wider( names_from = electrode, values_from = c( Sample, NMDS1, NMDS2 ) ) %>%
#|   
#|   filter( !is.na(NMDS1_H), !is.na(NMDS1_O) )
#| 
#| p_nmds <- ggplot(
#|   nmds_scores,
#|   aes( x = NMDS1, y = NMDS2, colour = treatment, shape = electrode )
#| ) +
#|   
#|   geom_segment( data = seg_nmds, aes( x = NMDS1_H, y = NMDS2_H, xend = NMDS1_O, yend = NMDS2_O ), inherit.aes = FALSE, linetype = 2, linewidth = .5 ) +
#|   
#|   geom_point( size = 4 ) +
#|   
#|   geom_text( aes( label = Sample ), vjust = -1, size = 3, show.legend = FALSE ) +
#|   
#|   geom_segment( data = vec_plot, aes( x = 0, y = 0, xend = NMDS1, yend = NMDS2 ), inherit.aes = FALSE, arrow = arrow( length = unit( .18, "cm" ) ), linewidth = .7 ) +
#|   
#|   geom_text( data = vec_plot, aes( x = NMDS1, y = NMDS2, label = Variable ), inherit.aes = FALSE, fontface = "bold", size = 3.5 ) +
#|   
#|   annotate(
#|     "text",
#|     x = Inf,
#|     y = Inf,
#|     hjust = 1.1,
#|     vjust = 1.3,
#|     label = paste0(
#|       "Stress = ",
#|       round(
#|         bac_nmds$stress,
#|         3
#|       )
#|     )
#|   ) +
#|   
#|   labs( x = "NMDS1", y = "NMDS2" ) +
#|   
#|   theme_classic( base_size = 11 )
#| 
#| # ============================================================
#| # 11. FIGURE 4B — UNIVARIATE dbRDA
#| # ============================================================
#| 
#| d_bray <- vegdist(
#|   bac_ra.v,
#|   method = "bray"
#| )
#| 
#| set.seed(123)
#| 
#| db_results <- map_dfr(
#|   env_vars,
#|   function(v) {
#|     
#|     form <- as.formula(
#|       paste0(
#|         "d_bray ~ ",
#|         v
#|       )
#|     )
#|     
#|     mod <- capscale(
#|       form,
#|       data = prok_sample
#|     )
#|     
#|     an <- anova(
#|       mod,
#|       permutations = 999
#|     )
#|     
#|     rsq <- RsquareAdj(
#|       mod
#|     )
#|     
#|     tibble(
#|       Variable = v,
#|       R2 = rsq$r.squared,
#|       Adj_R2 = rsq$adj.r.squared,
#|       F = an$F[1],
#|       p = an$`Pr(>F)`[1]
#|     )
#|   }
#| ) %>%
#|   
#|   mutate( FDR = p.adjust( p, method = "BH" ) )
#| 
#| cat("\n================ UNIVARIATE dbRDA ================\n")
#| 
#| print( db_results, n = Inf )
#| 
#| write.csv( db_results, "Table_univariate_dbRDA.csv", row.names = FALSE )
#| 
#| db_plot <- db_results %>%
#|   
#|   mutate(
#|     Variable_plot = recode(
#|       Variable,
#|       "N_mean" = "N",
#|       "C_mean" = "C",
#|       "H_mean" = "H",
#|       "S_mean" = "S",
#|       "MCA_mean" = "MCA",
#|       "MUF_mean" = "MUF"
#|     ),
#|     
#|     sig = case_when(
#|       FDR < .05 ~ "**",
#|       p < .05 ~ "*",
#|       TRUE ~ ""
#|     )
#|   )
#| 
#| p_dbrda <- ggplot(
#|   db_plot,
#|   aes( x = reorder( Variable_plot, R2 ), y = R2 )
#| ) +
#|   
#|   geom_col( width = .7 ) +
#|   
#|   geom_text( aes( label = sig ), hjust = -.3, size = 5 ) +
#|   
#|   coord_flip() +
#|   
#|   expand_limits(
#|     y = max(
#|       db_plot$R2
#|     ) * 1.15
#|   ) +
#|   
#|   labs( x = NULL, y = expression(R^2) ) +
#|   
#|   theme_classic( base_size = 11 )
#| 

# ===== ORIGINAL SOURCE: 08_Fig7_Mantel_Spearman.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 7411-7776; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Order environment Mantel/Spearman (original Fig.5 -> manuscript Fig.7)
#| 
#| # Combined Figure 4
#| p_environment <- p_nmds +
#|   p_dbrda +
#|   plot_layout(
#|     widths = c(
#|       1.3,
#|       1
#|     )
#|   )
#| 
#| p_environment
#| 
#| save_fig(
#|   p_environment,
#|   "Figure4_Community_environment",
#|   11,
#|   6
#| )
#| 
#| # ============================================================
#| # 12. ORDER TABLE FOR MANTEL ANALYSIS
#| # ============================================================
#| 
#| order_all <- order_abund %>%
#|   
#|   filter( Order != "Unclassified" ) %>%
#|   
#|   pivot_wider( names_from = Order, values_from = Abundance, values_fill = 0 ) %>%
#|   
#|   column_to_rownames(
#|     "Sample"
#|   ) %>%
#|   
#|   as.matrix()
#| 
#| # Exact sample alignment
#| order_all <- order_all[
#|   rownames(bac_ra.v),
#|   ,
#|   drop = FALSE
#| ]
#| 
#| # Filter:
#| # Order present in >=3 samples
#| # AND reaches >=1% in >=1 sample
#| keep_mantel <- apply(
#|   order_all,
#|   2,
#|   function(x) {
#|     
#|     sum(
#|       x > 0
#|     ) >= 3 &&
#|       max(x) >= 1
#|   }
#| )
#| 
#| order_tab <- order_all[
#|   ,
#|   keep_mantel,
#|   drop = FALSE
#| ]
#| 
#| cat("\nOrders included in Mantel analysis:\n")
#| 
#| print( colnames( order_tab ) )
#| 
#| # ============================================================
#| # 13. ORDER × ENVIRONMENT MANTEL + SPEARMAN
#| # ============================================================
#| 
#| env_mantel <- prok_sample[
#|   rownames(order_tab),
#|   env_vars,
#|   drop = FALSE
#| ]
#| 
#| set.seed(123)
#| 
#| taxa_env <- map_dfr(
#|   colnames(order_tab),
#|   function(ord) {
#|     
#|     map_dfr(
#|       env_vars,
#|       function(v) {
#|         
#|         x <- order_tab[
#|           ,
#|           ord
#|         ]
#|         
#|         e <- env_mantel[
#|           ,
#|           v
#|         ]
#|         
#|         # Mantel:
#|         # similarity of taxon-abundance distances
#|         # and environmental distances
#|         mt <- mantel(
#|           dist(x),
#|           dist(e),
#|           method = "spearman",
#|           permutations = 999
#|         )
#|         
#|         # Spearman:
#|         # direction of abundance-environment relationship
#|         sr <- suppressWarnings(
#|           cor.test(
#|             x,
#|             e,
#|             method = "spearman",
#|             exact = FALSE
#|           )
#|         )
#|         
#|         tibble(
#|           Order = ord,
#|           Variable = v,
#|           Mantel_r =
#|             unname(
#|               mt$statistic
#|             ),
#|           p =
#|             mt$signif,
#|           Spearman_rho =
#|             unname(
#|               sr$estimate
#|             )
#|         )
#|       }
#|     )
#|   }
#| )
#| 
#| # FDR within each environmental variable
#| taxa_env <- taxa_env %>%
#|   
#|   group_by(
#|     Variable
#|   ) %>%
#|   
#|   mutate( FDR_by_variable = p.adjust( p, method = "BH" ) ) %>%
#|   
#|   ungroup() %>%
#|   
#|   # Global correction across all Order × variable tests
#|   mutate( FDR_global = p.adjust( p, method = "BH" ) )
#| 
#| write.csv( taxa_env, "Table_Order_environment_Mantel.csv", row.names = FALSE )
#| 
#| cat("\n================ MANTEL p < 0.05 ================\n")
#| 
#| taxa_env %>%
#|   filter( p < .05 ) %>%
#|   arrange( p ) %>%
#|   print( n = Inf )
#| 
#| cat("\n================ WITHIN-VARIABLE FDR < 0.05 ================\n")
#| 
#| taxa_env %>%
#|   filter( FDR_by_variable < .05 ) %>%
#|   arrange( FDR_by_variable ) %>%
#|   print( n = Inf )
#| 
#| cat("\n================ GLOBAL FDR < 0.05 ================\n")
#| 
#| taxa_env %>%
#|   filter( FDR_global < .05 ) %>%
#|   arrange( FDR_global ) %>%
#|   print( n = Inf )
#| 
#| # ============================================================
#| # 14. FIGURE 5 — ORDER × ENVIRONMENT
#| # MANTEL r + SPEARMAN rho
#| # ============================================================
#| 
#| # Show Orders having at least one nominal Mantel p < 0.05
#| orders_keep_mantel <- taxa_env %>%
#|   filter( p < .05 ) %>%
#|   distinct( Order ) %>%
#|   pull(
#|     Order
#|   )
#| 
#| heat_data <- taxa_env %>%
#|   
#|   filter( Order %in% orders_keep_mantel ) %>%
#|   
#|   mutate(
#|     
#|     Variable = recode(
#|       as.character(
#|         Variable
#|       ),
#|       "N_mean" = "N",
#|       "C_mean" = "C",
#|       "H_mean" = "H",
#|       "S_mean" = "S",
#|       "MCA_mean" = "MCA",
#|       "MUF_mean" = "MUF"
#|     ),
#|     
#|     Variable = factor(
#|       Variable,
#|       levels = c(
#|         "pH",
#|         "ORP",
#|         "N",
#|         "C",
#|         "H",
#|         "S",
#|         "MCA",
#|         "MUF"
#|       )
#|     ),
#|     
#|     # Asterisk only for within-variable FDR < 0.05
#|     sig = ifelse(
#|       FDR_by_variable < .05,
#|       "*",
#|       ""
#|     )
#|   )
#| 
#| # Order taxa by environmental variable
#| # associated with their strongest Mantel r
#| order_levels_mantel <- heat_data %>%
#|   
#|   group_by(
#|     Order
#|   ) %>%
#|   
#|   slice_max(
#|     Mantel_r,
#|     n = 1,
#|     with_ties = FALSE
#|   ) %>%
#|   
#|   mutate( Variable_order = match( as.character( Variable ), c( "pH", "ORP", "N", "C", "H", "S", "MCA", "MUF" ) ) ) %>%
#|   
#|   arrange( desc( Variable_order ), Mantel_r ) %>%
#|   
#|   pull(
#|     Order
#|   )
#| 
#| heat_data$Order <- factor(
#|   heat_data$Order,
#|   levels = order_levels_mantel
#| )
#| 
#| p_mantel <- ggplot(
#|   heat_data,
#|   aes( x = Variable, y = Order )
#| ) +
#|   
#|   geom_tile( fill = "grey97", colour = "white", linewidth = .4 ) +
#|   
#|   geom_point( aes( size = Mantel_r, colour = Spearman_rho ) ) +
#|   
#|   geom_text( aes( label = sig ), colour = "black", fontface = "bold", size = 4 ) +
#|   
#|   scale_colour_gradient2( low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1,1), breaks = c( -1, -.5, 0, .5, 1 ), name = expression( "Spearman " * rho ) ) +
#|   
#|   scale_size_continuous( range = c( 1, 7 ), limits = c( 0, 1 ), breaks = c( 0, .25, .50, .75, 1 ), name = "Mantel r" ) +
#|   
#|   guides( size = guide_legend( order = 1 ), colour = guide_colourbar( order = 2 ) ) +
#|   
#|   labs( x = "Environmental variables", y = NULL ) +
#|   
#|   theme_classic( base_size = 11 ) +
#|   
#|   theme( axis.text.x = element_text( angle = 45, hjust = 1, face = "bold", colour = "black" ), axis.text.y = element_text( size = 8, colour = "black" ), legend.position = "right" )
#| 
#| p_mantel
#| 
#| save_fig(
#|   p_mantel,
#|   "Figure5_Order_environment",
#|   8.5,
#|   7.5
#| )
#| 
#| # ============================================================
#| # 15. FINAL SUMMARY
#| # ============================================================
#| 
#| cat("\n\n============================================================\n")
#| cat("ANALYSIS COMPLETED\n")
#| cat("============================================================\n")
#| 
#| cat("\nSamples analysed:\n")
#| print( sample_names( prok_data ) )
#| 
#| cat("\nNumber of samples:\n")
#| print( nsamples( prok_data ) )
#| 
#| cat("\nNumber of taxa/ASVs:\n")
#| print( ntaxa( prok_data ) )
#| 
#| cat("\nAlpha diversity:\n")
#| 
#| alpha_div %>%
#|   arrange(Sample) %>%
#|   as_tibble() %>%
#|   print(n = Inf) 
#| 
#| cat("\nStrong environmental correlations |rho| >= 0.60:\n")
#| cor_long %>%
#|   filter( abs(rho) >= .60 ) %>%
#|   print( n = Inf )
#| 
#| cat("\nUnivariate dbRDA:\n")
#| db_results %>%
#|   arrange( p ) %>%
#|   print( n = Inf )
#| 
#| cat("\nENVFIT:\n")
#| vec %>%
#|   arrange(p) %>%
#|   as_tibble() %>%
#|   print(n = Inf)
#| 
#| cat("\nV1HT1 Orders >=1%:\n")
#| HT1_taxa %>%
#|   filter( HT1 >= 1 ) %>%
#|   print( n = Inf )
#| 
#| cat("\nMantel associations — within-variable FDR <0.05:\n")
#| taxa_env %>%
#|   filter( FDR_by_variable < .05 ) %>%
#|   arrange( FDR_by_variable ) %>%
#|   print( n = Inf )
#| 
#| cat("\nGlobal Mantel FDR <0.05:\n")
#| taxa_env %>%
#|   filter( FDR_global < .05 ) %>%
#|   arrange( FDR_global ) %>%
#|   print( n = Inf )
#| 
#| cat("\nSession information:\n")
#| sessionInfo()
#| 
#| # ============================================================
#| # 16. FINAL EXPORT — INDIVIDUAL RESULT FIGURES
#| # ============================================================
#| # These exports preserve the six-result-figure structure used in the manuscript.
#| # Combined figures generated above are retained as optional alternatives.
#| 
#| save_fig(p_alpha,      "Figure1_Alpha_diversity",          8,   7)
#| save_fig(p_unifrac,    "Figure2_UniFrac_PCoA",            11,  5.5)
#| save_fig(p_order_bar,  "Figure3_Order_composition",        11,  6)
#| save_fig(p_order_heat, "Figure4_Order_Zscore_heatmap",     11,  7)
#| save_fig(p_environment,"Figure5_Community_environment",    11,  6)
#| save_fig(p_mantel,     "Figure6_Order_environment",        8.5, 7.5)
#| 
#| cat("\nFinal manuscript figure exports completed.\n")
#| ####################################################################################################
#| ###################################################################################################
#| # ============================================================
#| colnames(pathway)
#| ###############################################################################
#| ###############################################################################

# ===== ORIGINAL SOURCE: 09_Fig8_Tax4Fun_pathways.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 7778-8050; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Functional pathway plots (original Fig.7 -> manuscript Fig.8)
#| 
#| # FIGURE 7 - PREDICTED FUNCTIONAL POTENTIAL
#| # ============================================================
#| 
#| library(tidyverse)
#| library(readxl)
#| library(patchwork)
#| 
#| # ============================================================
#| # 1. IMPORT
#| # ============================================================
#| 
#| pathway <- read_excel("pathway_predition_all.xlsx")
#| functional <- read_excel("functional_predition_all.xlsx")
#| 
#| samples <- c("V1HT1","V1HT2","V1OT2","V2HT1","V2HT2","V2OT1","V2OT2","V3CT1","V3CT2")
#| 
#| sample_labels <- c(
#|   V1HT1="E-WI H T1", V1HT2="E-WI H T2", V1OT2="E-WI O T2",
#|   V2HT1="WI H T1",   V2HT2="WI H T2",
#|   V2OT1="WI O T1",   V2OT2="WI O T2",
#|   V3CT1="Control T1", V3CT2="Control T2"
#| )
#| 
#| cat("Pathways:", nrow(pathway), "\n")
#| cat("KOs:", nrow(functional), "\n")
#| 
#| stopifnot(all(samples %in% colnames(pathway)))
#| 
#| # ============================================================
#| # 2. CLEAN PATHWAY TABLE
#| # ============================================================
#| 
#| pathway <- pathway %>%
#|   mutate(
#|     level1 = as.character(level1),
#|     level2 = as.character(level2),
#|     level3 = as.character(level3),
#|     `Contrib. %` = as.numeric(`Contrib. %`),
#|     `Cumulative %` = as.numeric(`Cumulative %`)
#|   ) %>%
#|   filter(!is.na(level1), level1 != "")
#| 
#| # KEGG hierarchy:
#| # level1 = specific pathway
#| # level2 = functional subcategory
#| # level3 = broad functional category
#| 
#| 
#| # ============================================================
#| # 3. ROW Z-SCORE
#| # ============================================================
#| 
#| row_z <- function(x) {
#|   x <- as.numeric(x)
#|   s <- sd(x, na.rm = TRUE)
#|   if(is.na(s) || s == 0) return(rep(0, length(x)))
#|   (x - mean(x, na.rm = TRUE)) / s
#| }
#| 
#| 
#| # ============================================================
#| # 4. FIGURE 7A - TOP 20 SIMPER PATHWAYS
#| # ============================================================
#| 
#| top20 <- pathway %>%
#|   filter(!is.na(`Contrib. %`)) %>%
#|   arrange(desc(`Contrib. %`)) %>%
#|   slice_head(n = 20) %>%
#|   mutate(
#|     Pathway_label = paste0(
#|       level1, " (",
#|       sprintf("%.2f", `Contrib. %`), "%)"
#|     )
#|   )
#| 
#| pathway_order_A <- rev(top20$Pathway_label)
#| 
#| heat_A <- top20 %>%
#|   select(level1, Pathway_label, all_of(samples)) %>%
#|   pivot_longer(
#|     cols = all_of(samples),
#|     names_to = "Sample",
#|     values_to = "Abundance"
#|   ) %>%
#|   group_by(Pathway_label) %>%
#|   mutate(Zscore = row_z(Abundance)) %>%
#|   ungroup() %>%
#|   mutate(
#|     Sample = factor(
#|       Sample,
#|       levels = samples,
#|       labels = sample_labels[samples]
#|     ),
#|     Pathway_label = factor(
#|       Pathway_label,
#|       levels = pathway_order_A
#|     )
#|   )
#| 
#| p7A <- ggplot(
#|   heat_A,
#|   aes(x = Sample, y = Pathway_label, fill = Zscore)
#| ) +
#|   geom_tile(color = "white", linewidth = 0.15) +
#|   scale_fill_gradient2(
#|     low = "#2166AC",
#|     mid = "white",
#|     high = "#B2182B",
#|     midpoint = 0,
#|     name = "Z-score"
#|   ) +
#|   labs(
#|     title = "A. Pathways contributing most to predicted functional dissimilarity",
#|     x = NULL,
#|     y = NULL
#|   ) +
#|   theme_classic(base_size = 10) +
#|   theme(
#|     axis.text.x = element_text(angle = 45, hjust = 1),
#|     axis.text.y = element_text(size = 8),
#|     plot.title = element_text(face = "bold", size = 11),
#|     legend.position = "right"
#|   )
#| 
#| print(p7A)
#| 
#| 
#| # ============================================================
#| # 5. FIGURE 7B - TARGETED PATHWAYS
#| # ============================================================
#| 
#| target_pathways <- c(
#|   "ABC transporters",
#|   "Two-component system",
#|   "Quorum sensing",
#|   "Bacterial chemotaxis",
#|   "Fatty acid metabolism",
#|   "Fatty acid degradation",
#|   "Carbon metabolism",
#|   "Glycolysis / Gluconeogenesis",
#|   "Citrate cycle (TCA cycle)",
#|   "Butanoate metabolism",
#|   "Propanoate metabolism",
#|   "Benzoate degradation",
#|   "Degradation of aromatic compounds",
#|   "Oxidative phosphorylation",
#|   "Nitrogen metabolism",
#|   "Sulfur metabolism",
#|   "Methane metabolism"
#| )
#| 
#| found <- target_pathways[target_pathways %in% pathway$level1]
#| missing <- target_pathways[!target_pathways %in% pathway$level1]
#| 
#| cat("\nTarget pathways found:\n")
#| print(found)
#| 
#| cat("\nTarget pathways NOT found:\n")
#| print(missing)
#| 
#| target <- pathway %>%
#|   filter(level1 %in% target_pathways)
#| 
#| pathway_order_B <- rev(
#|   target_pathways[target_pathways %in% target$level1]
#| )
#| 
#| heat_B <- target %>%
#|   select(level1, all_of(samples)) %>%
#|   pivot_longer(
#|     cols = all_of(samples),
#|     names_to = "Sample",
#|     values_to = "Abundance"
#|   ) %>%
#|   group_by(level1) %>%
#|   mutate(Zscore = row_z(Abundance)) %>%
#|   ungroup() %>%
#|   mutate(
#|     Sample = factor(
#|       Sample,
#|       levels = samples,
#|       labels = sample_labels[samples]
#|     ),
#|     level1 = factor(
#|       level1,
#|       levels = pathway_order_B
#|     )
#|   )
#| 
#| p7B <- ggplot(
#|   heat_B,
#|   aes(x = Sample, y = level1, fill = Zscore)
#| ) +
#|   geom_tile(color = "white", linewidth = 0.15) +
#|   scale_fill_gradient2(
#|     low = "#2166AC",
#|     mid = "white",
#|     high = "#B2182B",
#|     midpoint = 0,
#|     name = "Z-score"
#|   ) +
#|   labs(
#|     title = "B. Predicted pathways relevant to the experimental framework",
#|     x = NULL,
#|     y = NULL
#|   ) +
#|   theme_classic(base_size = 10) +
#|   theme(
#|     axis.text.x = element_text(angle = 45, hjust = 1),
#|     axis.text.y = element_text(size = 8),
#|     plot.title = element_text(face = "bold", size = 11),
#|     legend.position = "right"
#|   )
#| 
#| print(p7B)
#| 
#| 
#| # ============================================================
#| # 6. COMBINED FIGURE 7
#| # ============================================================
#| 
#| Figure7 <- p7A / p7B +
#|   plot_layout(heights=c(1.05, 1))
#| 
#| print(Figure7)
#| 
#| 
#| # ============================================================
#| # 7. EXPORT FIGURES
#| # ============================================================
#| 
#| ggsave("Figure7A_Functional_SIMPER.pdf", p7A, width=11, height=8)
#| ggsave("Figure7A_Functional_SIMPER.svg", p7A, width=11, height=8)
#| 
#| ggsave("Figure7B_Functional_Targeted.pdf", p7B, width=11, height=7)
#| ggsave("Figure7B_Functional_Targeted.svg", p7B, width=11, height=7)
#| 
#| ggsave("Figure7_Functional_prediction.pdf", Figure7, width=11, height=15)
#| ggsave("Figure7_Functional_prediction.svg", Figure7, width=11, height=15)
#| 
#| 
#| # ============================================================
#| # 8. EXPORT SOURCE TABLES
#| # ============================================================
#| 
#| write.csv(
#|   top20,
#|   "Table_Figure7A_TOP20_SIMPER.csv",
#|   row.names=FALSE
#| )
#| 
#| write.csv(
#|   target,
#|   "Table_Figure7B_Target_pathways.csv",
#|   row.names=FALSE
#| )
#| 
#| 
#| # ============================================================
#| # 9. SUMMARY
#| # ============================================================
#| 
#| cat("\n========================================\n")
#| cat("FIGURE 7 SUMMARY\n")
#| cat("========================================\n")
#| cat("Total pathways:", nrow(pathway), "\n")
#| cat("Total KOs:", nrow(functional), "\n")
#| cat("Figure 7A pathways:", nrow(top20), "\n")
#| cat("Figure 7B pathways:", nrow(target), "\n")
#| cat("========================================\n")
#| 
#| head(
#|   pathway %>%

# ===== ORIGINAL SOURCE: 11_Analisi_ambientali_seconda_versione.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 16971-17830; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Later original environment ordination and dbRDA exploration
#| 
#| # END FIGURE 3
#| ###############################################################################
#| ###############################################################################
#| # END
#| ###############################################################################
#| ###############################################################################
#| ###############################################################################
#| #########################################env fit con le variabili ambientali scalate################################
#| library(phyloseq)
#| library(vegan)
#| library(dplyr)
#| # Trasforma in percentuale (relative abundance)
#| bac_ra <- transform_sample_counts(prok_data, function(x) { (x / sum(x)) * 100 })
#| 
#| # Funzione per convertire OTU table in matrice compatibile vegan
#| psotu2veg <- function(physeq) {
#|   OTU <- otu_table(physeq)
#|   if (taxa_are_rows(OTU)) OTU <- t(OTU)
#|   return(as(OTU, "matrix"))
#| }
#| 
#| # OTU table per vegan
#| bac_ra.v <- psotu2veg(bac_ra)
#| #2️⃣ nMDS con distanze Jaccard
#| 
#| 
#| set.seed(123)  # per riproducibilità
#| bac.v_nmds_bray <- metaMDS(
#|   bac_ra.v,
#|   distance = "bray",
#|   trymax = 50
#| )
#| 
#| # ============================================================
#| # 3. PREPARAZIONE DATI AMBIENTALI
#| # ============================================================
#| 
#| env_data <- data.frame(sample_data(prok_data))[, c(
#|   "Samples",
#|   "treatment",
#|   "area",
#|   "time",
#|   "pH",
#|   "ORP",
#|   "N_mean",
#|   "C_mean",
#|   "H_mean",
#|   "S_mean",
#|   "MCA_mean",
#|   "MUF_mean"
#| )]
#| 
#| # Controllo iniziale
#| str(env_data)
#| 
#| 
#| # ============================================================
#| # 4. CONVERSIONE DELLE VARIABILI AMBIENTALI IN NUMERICO
#| # ============================================================
#| 
#| # Queste colonne sono state importate come character
#| # perché nel file i decimali sono indicati con la virgola
#| 
#| env_vars <- c(
#|   "pH",
#|   "ORP",
#|   "N_mean",
#|   "C_mean",
#|   "H_mean",
#|   "S_mean",
#|   "MCA_mean",
#|   "MUF_mean"
#| )
#| 
#| env_data[env_vars] <- lapply(
#|   env_data[env_vars],
#|   function(x) as.numeric(gsub(",", ".", x))
#| )
#| 
#| # Adesso controlliamo nuovamente
#| str(env_data)
#| 
#| 
#| # ============================================================
#| # 5. SELEZIONE DELLE VARIABILI PER ENVFIT
#| # ============================================================
#| 
#| env_data_num <- env_data[, env_vars]
#| 
#| colnames(env_data_num)
#| 
#| # Controllo eventuali NA
#| colSums(is.na(env_data_num))
#| 
#| 
#| # ============================================================
#| # 6. CONTROLLO ORDINE DEI CAMPIONI
#| # ============================================================
#| 
#| all(
#|   rownames(env_data_num) ==
#|     rownames(bac_ra.v)
#| )
#| 
#| # Se restituisce FALSE, riordina:
#| if (!all(rownames(env_data_num) == rownames(bac_ra.v))) {
#|   
#|   env_data_num <- env_data_num[
#|     rownames(bac_ra.v),
#|     ,
#|     drop = FALSE
#|   ]
#| }
#| 
#| 
#| # ============================================================
#| # 7. SCALING DELLE VARIABILI AMBIENTALI
#| # ============================================================
#| 
#| env_scaled <- as.data.frame(
#|   scale(env_data_num)
#| )
#| 
#| rownames(env_scaled) <- rownames(env_data_num)
#| 
#| # Controllo
#| head(env_scaled)
#| 
#| 
#| # ============================================================
#| # 8. ENVFIT SULLA nMDS BRAY-CURTIS
#| # ============================================================
#| 
#| set.seed(123)
#| 
#| env_fitting_scaled <- envfit(
#|   bac.v_nmds_bray,
#|   env_scaled,
#|   permutations = 9999
#| )
#| 
#| # Risultati
#| env_fitting_scaled
#| 
#| ***VECTORS
#| 
#| NMDS1    NMDS2     r2 Pr(>r)   
#| pH        0.34532 -0.93848 0.3822 0.2198   
#| ORP      -0.63613 -0.77158 0.9214 0.0011 **
#|   N_mean    0.78097  0.62457 0.1744 0.5593   
#| C_mean   -0.70357 -0.71063 0.6921 0.0297 * 
#|   H_mean   -0.16147  0.98688 0.2052 0.5102   
#| S_mean   -0.20112  0.97957 0.4687 0.1377   
#| MCA_mean -0.32506 -0.94569 0.2756 0.3279   
#| MUF_mean  0.12407 -0.99227 0.3115 0.3157   
#| ---
#|   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#| Permutation: free
#| Number of permutations: 9999
#| 
#| 
#| 
#| 
#| # Visualizza solo vettori significativi p < 0.05
#| sig_vars <- names(env_fitting_scaled$vectors$pvals)[env_fitting_scaled$vectors$pvals < 0.05]
#| env_fitting_sig <- envfit(bac.v_nmds_bray, env_scaled[, sig_vars], perm = 9999)
#| #⃣ Plot nMDS con vettori ambientali scalati
#| # Plot base nMDS
#| plot(bac.v_nmds_j, display = "sites", main = "nMDS bray - Env Fit Scaled")
#| 
#| # Aggiungi vettori significativi
#| plot(env_fitting_sig, col = "blue", p.max = 0.05) 
#| 
#| plot(bac.v_nmds_bray, display = "sites", main = "nMDS Jaccard_Cl")
#| text(bac.v_nmds_bray, display = "sites",
#|      labels = sample_names(prok_data),
#|      cex = 0.7, pos = 3, offset = 0.4)
#| with(data.frame(sample_data(prok_data)), ordisurf(bac.v_nmds_bray, ORP, add = TRUE, col = "green4"))
#| plot(env_fitting, p.max = 0.1, col = "red") 
#| ######################################################################################
#| # ============================================================
#| # nMDS BRAY-CURTIS + ENVFIT + ORDISURF
#| # ============================================================
#| 
#| library(phyloseq)
#| library(vegan)
#| library(dplyr)
#| 
#| # ============================================================
#| # 1. ABBONDANZA RELATIVA
#| # ============================================================
#| 
#| bac_ra <- transform_sample_counts(
#|   prok_data,
#|   function(x) {
#|     (x / sum(x)) * 100
#|   }
#| )
#| 
#| 
#| # ============================================================
#| # 2. CONVERSIONE PHYLOSEQ -> MATRICE VEGAN
#| # ============================================================
#| 
#| psotu2veg <- function(physeq) {
#|   
#|   OTU <- otu_table(physeq)
#|   
#|   if (taxa_are_rows(OTU)) {
#|     OTU <- t(OTU)
#|   }
#|   
#|   as(OTU, "matrix")
#| }
#| 
#| bac_ra.v <- psotu2veg(bac_ra)
#| 
#| 
#| # ============================================================
#| # 3. nMDS BRAY-CURTIS
#| # ============================================================
#| 
#| set.seed(123)
#| 
#| bac.v_nmds_bray <- metaMDS(
#|   bac_ra.v,
#|   distance = "bray",
#|   trymax = 100,
#|   autotransform = FALSE
#| )
#| 
#| # Stress della nMDS
#| bac.v_nmds_bray$stress
#| 
#| 
#| # ============================================================
#| # 4. PREPARAZIONE METADATA AMBIENTALI
#| # ============================================================
#| 
#| env_data <- data.frame(
#|   sample_data(prok_data)
#| )[, c(
#|   "Samples",
#|   "treatment",
#|   "area",
#|   "time",
#|   "pH",
#|   "ORP",
#|   "N_mean",
#|   "C_mean",
#|   "H_mean",
#|   "S_mean",
#|   "MCA_mean",
#|   "MUF_mean"
#| )]
#| 
#| # Controllo iniziale
#| str(env_data)
#| 
#| 
#| # ============================================================
#| # 5. CONVERSIONE DELLE VARIABILI AMBIENTALI IN NUMERICO
#| # ============================================================
#| 
#| env_vars <- c(
#|   "pH",
#|   "ORP",
#|   "N_mean",
#|   "C_mean",
#|   "H_mean",
#|   "S_mean",
#|   "MCA_mean",
#|   "MUF_mean"
#| )
#| 
#| env_data[env_vars] <- lapply(
#|   env_data[env_vars],
#|   function(x) {
#|     as.numeric(
#|       gsub(",", ".", x)
#|     )
#|   }
#| )
#| 
#| # Controllo
#| str(env_data)
#| 
#| 
#| # ============================================================
#| # 6. DATAFRAME SOLO CON VARIABILI NUMERICHE
#| # ============================================================
#| 
#| env_data_num <- env_data[
#|   ,
#|   env_vars,
#|   drop = FALSE
#| ]
#| 
#| 
#| # ============================================================
#| # 7. CONTROLLO EVENTUALI NA
#| # ============================================================
#| 
#| cat("\nNumero di NA per variabile:\n")
#| 
#| print(
#|   colSums(
#|     is.na(env_data_num)
#|   )
#| )
#| 
#| 
#| # ============================================================
#| # 8. SINCRONIZZAZIONE ORDINE CAMPIONI
#| # ============================================================
#| 
#| cat("\nControllo ordine campioni:\n")
#| 
#| print(
#|   all(
#|     rownames(env_data_num) ==
#|       rownames(bac_ra.v)
#|   )
#| )
#| 
#| # Riordina comunque i metadata nello stesso ordine della OTU table
#| env_data_num <- env_data_num[
#|   rownames(bac_ra.v),
#|   ,
#|   drop = FALSE
#| ]
#| 
#| env_data <- env_data[
#|   rownames(bac_ra.v),
#|   ,
#|   drop = FALSE
#| ]
#| 
#| 
#| # ============================================================
#| # 9. ELIMINA EVENTUALI VARIABILI A VARIANZA ZERO
#| # ============================================================
#| 
#| sd_vars <- apply(
#|   env_data_num,
#|   2,
#|   sd,
#|   na.rm = TRUE
#| )
#| 
#| cat("\nDeviazione standard variabili ambientali:\n")
#| 
#| print(sd_vars)
#| 
#| env_data_num <- env_data_num[
#|   ,
#|   sd_vars > 0,
#|   drop = FALSE
#| ]
#| 
#| 
#| # ============================================================
#| # 10. SCALING DELLE VARIABILI AMBIENTALI
#| # ============================================================
#| 
#| env_scaled <- as.data.frame(
#|   scale(env_data_num)
#| )
#| 
#| rownames(env_scaled) <- rownames(env_data_num)
#| 
#| 
#| # ============================================================
#| # 11. ENVFIT
#| # ============================================================
#| 
#| set.seed(123)
#| 
#| env_fitting_scaled <- envfit(
#|   bac.v_nmds_bray,
#|   env_scaled,
#|   permutations = 9999
#| )
#| 
#| cat("\nRISULTATI ENVFIT:\n")
#| 
#| print(env_fitting_scaled)
#| 
#| 
#| # ============================================================
#| # 12. TABELLA RISULTATI ENVFIT
#| # ============================================================
#| 
#| envfit_results <- data.frame(
#|   Variable = rownames(
#|     env_fitting_scaled$vectors$arrows
#|   ),
#|   
#|   R2 = env_fitting_scaled$vectors$r,
#|   
#|   P_value = env_fitting_scaled$vectors$pvals
#| )
#| 
#| envfit_results <- envfit_results %>%
#|   mutate(
#|     Significance = case_when(
#|       P_value < 0.001 ~ "***",
#|       P_value < 0.01  ~ "**",
#|       P_value < 0.05  ~ "*",
#|       P_value < 0.10  ~ ".",
#|       TRUE            ~ "ns"
#|     )
#|   ) %>%
#|   arrange(P_value)
#| 
#| cat("\nTABELLA ENVFIT:\n")
#| 
#| print(envfit_results)
#| 
#| 
#| # ============================================================
#| # 13. SELEZIONE DELLE VARIABILI PIU' SIGNIFICATIVE
#| # ============================================================
#| 
#| # Soglia esplorativa
#| p_threshold <- 0.10
#| 
#| sig_vars <- envfit_results %>%
#|   filter(P_value < p_threshold) %>%
#|   pull(Variable)
#| 
#| cat("\nVariabili con p <", p_threshold, ":\n")
#| 
#| print(sig_vars)
#| 
#| 
#| # ============================================================
#| # 14. GRAFICO nMDS BRAY-CURTIS
#| # ============================================================
#| 
#| plot(
#|   bac.v_nmds_bray,
#|   display = "sites",
#|   type = "n",
#|   main = paste0(
#|     "nMDS Bray-Curtis | Stress = ",
#|     round(
#|       bac.v_nmds_bray$stress,
#|       3
#|     )
#|   )
#| )
#| 
#| points(
#|   bac.v_nmds_bray,
#|   display = "sites",
#|   pch = 19,
#|   cex = 1.5
#| )
#| 
#| text(
#|   bac.v_nmds_bray,
#|   display = "sites",
#|   labels = sample_names(prok_data),
#|   cex = 0.75,
#|   pos = 3,
#|   offset = 0.4
#| )
#| 
#| 
#| # ============================================================
#| # 15. AGGIUNTA DEI VETTORI ENVFIT
#| # ============================================================
#| 
#| plot(
#|   env_fitting_scaled,
#|   p.max = p_threshold,
#|   col = "red"
#| )
#| 
#| 
#| # ============================================================
#| # 16. ORDISURF PER LE VARIABILI PIU' SIGNIFICATIVE
#| # ============================================================
#| 
#| ordisurf_results <- list()
#| 
#| if (length(sig_vars) > 0) {
#|   
#|   for (var in sig_vars) {
#|     
#|     cat(
#|       "\n====================================\n",
#|       "ORDISURF:", var,
#|       "\n====================================\n"
#|     )
#|     
#|     # Grafico vuoto
#|     plot(
#|       bac.v_nmds_bray,
#|       display = "sites",
#|       type = "n",
#|       main = paste(
#|         "nMDS Bray-Curtis -",
#|         var
#|       )
#|     )
#|     
#|     # Campioni
#|     points(
#|       bac.v_nmds_bray,
#|       display = "sites",
#|       pch = 19,
#|       cex = 1.5
#|     )
#|     
#|     # Etichette campioni
#|     text(
#|       bac.v_nmds_bray,
#|       display = "sites",
#|       labels = sample_names(prok_data),
#|       cex = 0.75,
#|       pos = 3,
#|       offset = 0.4
#|     )
#|     
#|     # Superficie ambientale
#|     fit <- ordisurf(
#|       bac.v_nmds_bray,
#|       env_data_num[[var]],
#|       add = TRUE,
#|       knots = 3
#|     )
#|     
#|     # Salva il modello
#|     ordisurf_results[[var]] <- fit
#|     
#|     print(fit)
#|   }
#|   
#| } else {
#|   
#|   cat(
#|     "\nNessuna variabile con p <",
#|     p_threshold,
#|     "\n"
#|   )
#| }
#| 
#| 
#| # ============================================================
#| # 17. SUMMARY DEI MODELLI ORDISURF
#| # ============================================================
#| 
#| if (length(ordisurf_results) > 0) {
#|   
#|   cat("\nSUMMARY ORDISURF:\n")
#|   
#|   lapply(
#|     ordisurf_results,
#|     summary
#|   )
#| }
#| ################################################################################
#| # ============================================================
#| # dbRDA / CAP SU BRAY-CURTIS
#| # ============================================================
#| 
#| library(vegan)
#| 
#| # Matrice comunità
#| comm <- bac_ra.v
#| 
#| # Variabili ambientali numeriche già convertite
#| env <- env_data_num
#| 
#| # Controllo ordine campioni
#| env <- env[rownames(comm), , drop = FALSE]
#| 
#| # Modello completo
#| dbrda_full <- capscale(
#|   comm ~ pH + ORP + N_mean + C_mean + H_mean + S_mean +
#|     MCA_mean + MUF_mean,
#|   data = env,
#|   distance = "bray",
#|   add = TRUE
#| )
#| 
#| # Test globale
#| anova(dbrda_full, permutations = 9999)
#| 
#| # Test per singola variabile
#| anova(dbrda_full, by = "terms", permutations = 9999)
#| 
#| # Test degli assi canonici
#| anova(dbrda_full, by = "axis", permutations = 9999)
#| 
#| # Percentuale di varianza spiegata
#| RsquareAdj(dbrda_full)
#| 
#| # Grafico
#| plot(
#|   dbrda_full,
#|   display = c("sites", "bp"),
#|   main = "dbRDA - Bray-Curtis"
#| )
#| 
#| text(
#|   dbrda_full,
#|   display = "sites",
#|   labels = rownames(comm),
#|   cex = 0.8,
#|   pos = 3
#| )
#| ################################################################################
#| cor_env <- cor(
#|   env_data_num,
#|   use = "pairwise.complete.obs",
#|   method = "spearman"
#| )
#| 
#| round(cor_env, 2)
#| > round(cor_env, 2)
#| pH   ORP N_mean C_mean H_mean S_mean MCA_mean MUF_mean
#| pH        1.00  0.27  -0.64   0.03  -0.38  -0.22     0.69     0.13
#| ORP       0.27  1.00   0.10   0.88   0.07  -0.25     0.72     0.52
#| N_mean   -0.64  0.10   1.00   0.38   0.65  -0.10    -0.17     0.22
#| C_mean    0.03  0.88   0.38   1.00   0.22  -0.25     0.48     0.47
#| H_mean   -0.38  0.07   0.65   0.22   1.00   0.57    -0.30    -0.23
#| S_mean   -0.22 -0.25  -0.10  -0.25   0.57   1.00    -0.63    -0.77
#| MCA_mean  0.69  0.72  -0.17   0.48  -0.30  -0.63     1.00     0.67
#| MUF_mean  0.13  0.52   0.22   0.47  -0.23  -0.77     0.67     1.00
#| 
#| vif.cca(dbrda_full)
#| 
#| save.image()
#| 
#| envfit_results
#| 
#| 
#| > envfit_results
#| Variable         R2 P_value Significance
#| ORP           ORP 0.97782217  0.0007          ***
#|   C_mean     C_mean 0.73970009  0.0235            *
#|   S_mean     S_mean 0.49217639  0.1175           ns
#| MUF_mean MUF_mean 0.46549601  0.1450           ns
#| MCA_mean MCA_mean 0.39277531  0.1616           ns
#| pH             pH 0.38390946  0.2097           ns
#| H_mean     H_mean 0.15052680  0.6261           ns
#| N_mean     N_mean 0.07967847  0.8206           ns
#| 
#| 
#| # ============================================================
#| # MODELLO 1: CONDIZIONI FISICO-CHIMICHE
#| # ============================================================
#| 
#| dbrda_env <- capscale(
#|   comm ~ pH + ORP,
#|   data = env,
#|   distance = "bray",
#|   add = TRUE
#| )
#| 
#| vif.cca(dbrda_env)
#| 
#| anova(
#|   dbrda_env,
#|   permutations = 9999
#| )
#| 
#| anova(
#|   dbrda_env,
#|   by = "margin",
#|   permutations = 9999
#| )
#| 
#| RsquareAdj(dbrda_env)
#| 
#| 
#|  # MODELLO 1: CONDIZIONI FISICO-CHIMICHE
#| > # ============================================================
#| > 
#| > dbrda_env <- capscale(
#| +   comm ~ pH + ORP,
#| +   data = env,
#| +   distance = "bray",
#| +   add = TRUE
#| + )
#| > 
#| > vif.cca(dbrda_env)
#|       pH      ORP 
#| 1.003684 1.003684 
#| > 
#| > anova(
#| +   dbrda_env,
#| +   permutations = 9999
#| + )
#| Permutation test for capscale under reduced model
#| Permutation: free
#| Number of permutations: 9999
#| 
#| Model: capscale(formula = comm ~ pH + ORP, data = env, distance = "bray", add = TRUE)
#|          Df SumOfSqs      F Pr(>F)    
#| Model     2   1.1458 1.3696  3e-04 ***
#| Residual  6   2.5097                  
#| ---
#| Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#| > 
#| > anova(
#| +   dbrda_env,
#| +   by = "margin",
#| +   permutations = 9999
#| + )
#| Permutation test for capscale under reduced model
#| Marginal effects of terms
#| Permutation: free
#| Number of permutations: 9999
#| 
#| Model: capscale(formula = comm ~ pH + ORP, data = env, distance = "bray", add = TRUE)
#|          Df SumOfSqs      F Pr(>F)    
#| pH        1  0.52130 1.2463 0.0169 *  
#| ORP       1  0.61793 1.4773 0.0002 ***
#| Residual  6  2.50969                  
#| ---
#| Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#| > 
#| > RsquareAdj(dbrda_env)
#| $r.squared
#| [1] 0.3134435
#| 
#| $adj.r.squared
#| [1] 0.08459136
#| 
#| 
#| save.image()
#| 
#| 
#| # ============================================================
#| # MODELLO 2: COMPOSIZIONE DEL SEDIMENTO
#| # ============================================================
#| 
#| dbrda_element <- capscale(
#|   comm ~ C_mean + N_mean,
#|   data = env,
#|   distance = "bray",
#|   add = TRUE
#| )
#| 
#| vif.cca(dbrda_element)
#| 
#| anova(
#|   dbrda_element,
#|   permutations = 9999
#| )
#| 
#| anova(
#|   dbrda_element,
#|   by = "margin",
#|   permutations = 9999
#| )
#| 
#| RsquareAdj(dbrda_element)
#| 
#| ################################################################################
#| 
#| save.image()
#| # ============================================================
#| # MODELLO 3: ATTIVITA' ENZIMATICHE
#| # ============================================================
#| 
#| dbrda_enz <- capscale(
#|   comm ~ MCA_mean + MUF_mean,
#|   data = env,
#|   distance = "bray",
#|   add = TRUE
#| )
#| 
#| vif.cca(dbrda_enz)
#| 
#| anova(
#|   dbrda_enz,
#|   permutations = 9999
#| )
#| 
#| anova(
#|   dbrda_enz,
#|   by = "margin",
#|   permutations = 9999
#| )
#| 
#| RsquareAdj(dbrda_enz)
#| #################################################################################
#| dbrda_ORP <- capscale(
#|   comm ~ ORP,
#|   data = env,
#|   distance = "bray",
#|   add = TRUE
#| )
#| 
#| dbrda_C <- capscale(
#|   comm ~ C_mean,
#|   data = env,
#|   distance = "bray",
#|   add = TRUE
#| )
#| 
#| anova(dbrda_ORP, permutations = 9999)
#| anova(dbrda_C, permutations = 9999)
#| 
#| RsquareAdj(dbrda_ORP)
#| RsquareAdj(dbrda_C)
#| #################################################################################
#| 
#| anova(dbrda_env, by = "margin", permutations = 9999)
#| RsquareAdj(dbrda_env)
#| 
#| anova(dbrda_element, by = "margin", permutations = 9999)
#| RsquareAdj(dbrda_element)
#| 
#| anova(dbrda_enz, by = "margin", permutations = 9999)
#| RsquareAdj(dbrda_enz)
#| ####
#| > anova(dbrda_enz, by = "margin", permutations = 9999)
#| Permutation test for capscale under reduced model
#| Marginal effects of terms
#| Permutation: free
#| Number of permutations: 9999
#| 
#| Model: capscale(formula = comm ~ MCA_mean + MUF_mean, data = env, distance = "bray", add = TRUE)
#| Df SumOfSqs      F Pr(>F)  
#| MCA_mean  1  0.56268 1.3224 0.0156 *
#|   MUF_mean  1  0.57756 1.3574 0.0113 *
#|   Residual  6  2.55290                
#| ---
#|   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#| > RsquareAdj(dbrda_enz)
#| $r.squared
#| [1] 0.3016223
#| 
#| $adj.r.squared
#| [1] 0.06882976
#| 
#| ###############################################################################
#| # ============================================================
#| # MANTEL TEST
#| # TAXA vs VARIABILI AMBIENTALI
#| # ============================================================
#| 
#| library(phyloseq)
#| library(vegan)
#| library(dplyr)
#| library(tidyr)
#| 
#| 
#| # ============================================================
#| # 1. ABBONDANZA RELATIVA
#| # ============================================================
#| 
#| prok_ra <- transform_sample_counts(
#|   prok_data,

# ===== ORIGINAL SOURCE: 12_Mantel_finale_seconda_versione.R =====
#| # ARCHIVAL EXTRACTION FROM ITEM4_24_26_(xarticolo).R
#| # Original source lines 18390-19377; unchanged below.
#| # NOTE: archival code may depend on objects created in other sections and may contain console output.
#| # Later original Mantel/Spearman plotting exploration
#| 
#| ############################################################
#| 
#| colSums(
#|   is.na(env)
#| )
#| 
#| 
#| ############################################################
#| # Rimuoviamo eventuali campioni con NA ambientali
#| ############################################################
#| 
#| keep_rows <- complete.cases(env)
#| 
#| df_sel <- df_sel[
#|   keep_rows,
#|   ,
#|   drop = FALSE
#| ]
#| 
#| env <- env[
#|   keep_rows,
#|   ,
#|   drop = FALSE
#| ]
#| 
#| 
#| ############################################################
#| # 6. RIMOZIONE VARIABILI CON SD = 0
#| ############################################################
#| 
#| sd_env <- apply(
#|   env,
#|   2,
#|   sd,
#|   na.rm = TRUE
#| )
#| 
#| 
#| print(sd_env)
#| 
#| 
#| keep_env <- names(
#|   sd_env[
#|     is.finite(sd_env) &
#|       sd_env > 0
#|   ]
#| )
#| 
#| 
#| env <- env[
#|   ,
#|   keep_env,
#|   drop = FALSE
#| ]
#| 
#| 
#| ############################################################
#| # 7. SCALING VARIABILI AMBIENTALI
#| ############################################################
#| 
#| env_scaled <- as.data.frame(
#|   scale(env)
#| )
#| 
#| 
#| rownames(env_scaled) <- rownames(env)
#| 
#| 
#| ############################################################
#| # CONTROLLO
#| ############################################################
#| 
#| summary(env_scaled)
#| 
#| 
#| ############################################################
#| # 8. CONTROLLO ORDINE CAMPIONI
#| ############################################################
#| 
#| stopifnot(
#|   identical(
#|     rownames(df_sel),
#|     rownames(env_scaled)
#|   )
#| )
#| 
#| 
#| ############################################################
#| # 9. SELEZIONE TAXA PER Mantel_test
#| ############################################################
#| #
#| # Ogni Order viene analizzato separatamente.
#| #
#| ############################################################
#| 
#| taxa_names <- colnames(df_sel)
#| 
#| 
#| spec_select_list <- setNames(
#|   as.list(
#|     seq_along(taxa_names)
#|   ),
#|   taxa_names
#| )
#| 
#| 
#| ############################################################
#| # 10. MANTEL TEST
#| ############################################################
#| #
#| # Distanze:
#| # taxa univariati -> Euclidean
#| # variabili ambientali scalate -> Euclidean
#| #
#| ############################################################
#| 
#| set.seed(123)
#| 
#| 
#| mantel <- linkET::mantel_test(
#|   
#|   spec = df_sel,
#|   
#|   env = env_scaled,
#|   
#|   spec_dist = "euclidean",
#|   
#|   env_dist = "euclidean",
#|   
#|   spec_select = spec_select_list
#|   
#| ) %>%
#|   
#|   mutate(
#|     
#|     ##########################################################
#|     # CLASSI Mantel r
#|     ##########################################################
#|     
#|     rd = cut(
#|       r,
#|       breaks = c(
#|         -Inf,
#|         0.2,
#|         0.4,
#|         Inf
#|       ),
#|       labels = c(
#|         "< 0.2",
#|         "0.2 - 0.4",
#|         ">= 0.4"
#|       )
#|     ),
#|     
#|     ##########################################################
#|     # CLASSI p-value
#|     ##########################################################
#|     
#|     pd = cut(
#|       p,
#|       breaks = c(
#|         -Inf,
#|         0.01,
#|         0.05,
#|         Inf
#|       ),
#|       labels = c(
#|         "< 0.01",
#|         "0.01 - 0.05",
#|         ">= 0.05"
#|       )
#|     )
#|   )
#| 
#| 
#| ############################################################
#| # 11. FDR / BENJAMINI-HOCHBERG
#| ############################################################
#| #
#| # La conserviamo nella tabella.
#| #
#| # Per il grafico, per ora, usiamo il p nominale,
#| # come nel grafico DJP.
#| #
#| ############################################################
#| 
#| mantel <- mantel %>%
#|   
#|   group_by(env) %>%
#|   
#|   mutate(
#|     FDR = p.adjust(
#|       p,
#|       method = "BH"
#|     )
#|   ) %>%
#|   
#|   ungroup()
#| 
#| 
#| ############################################################
#| # 12. ORDINE RISULTATI
#| ############################################################
#| 
#| mantel <- mantel %>%
#|   
#|   arrange(
#|     p,
#|     desc(r)
#|   )
#| 
#| 
#| ############################################################
#| # 13. VISUALIZZAZIONE RISULTATI
#| ############################################################
#| 
#| print(mantel)
#| 
#| 
#| ############################################################
#| # 14. SALVATAGGIO RISULTATI COMPLETI
#| ############################################################
#| 
#| write_xlsx(
#|   mantel,
#|   "Mantel_Order_environment_ALL_results.xlsx"
#| )
#| 
#| 
#| ############################################################
#| # 15. TAXA CON ALMENO UNA ASSOCIAZIONE SIGNIFICATIVA
#| ############################################################
#| 
#| taxa_sig <- mantel %>%
#|   
#|   filter(
#|     p < 0.05
#|   ) %>%
#|   
#|   pull(spec) %>%
#|   
#|   unique()
#| 
#| 
#| cat(
#|   "\nOrder con almeno una associazione p < 0.05:\n"
#| )
#| 
#| print(taxa_sig)
#| 
#| 
#| ############################################################
#| # 16. RISULTATI SIGNIFICATIVI
#| ############################################################
#| 
#| mantel_sig <- mantel %>%
#|   
#|   filter(
#|     spec %in% taxa_sig,
#|     p < 0.05
#|   )
#| 
#| 
#| ############################################################
#| # 17. RISULTATI NON SIGNIFICATIVI
#| #
#| # SOLO per i taxa che possiedono almeno una relazione
#| # significativa.
#| #
#| ############################################################
#| 
#| mantel_nonsig <- mantel %>%
#|   
#|   filter(
#|     spec %in% taxa_sig,
#|     p >= 0.05
#|   )
#| 
#| 
#| ############################################################
#| # 18. TABELLA SOLO SIGNIFICATIVI
#| ############################################################
#| 
#| mantel_significant_table <- mantel %>%
#|   
#|   filter(
#|     p < 0.05
#|   ) %>%
#|   
#|   arrange(
#|     p,
#|     desc(r)
#|   )
#| 
#| 
#| print(
#|   mantel_significant_table
#| )
#| 
#| 
#| ############################################################
#| # 19. SALVA TABELLA SIGNIFICATIVI
#| ############################################################
#| 
#| write_xlsx(
#|   mantel_significant_table,
#|   "Mantel_Order_environment_significant.xlsx"
#| )
#| 
#| 
#| ############################################################
#| # 20. RISULTATI CHE SUPERANO FDR
#| ############################################################
#| 
#| mantel_FDR <- mantel %>%
#|   
#|   filter(
#|     FDR < 0.05
#|   )
#| 
#| 
#| print(
#|   mantel_FDR
#| )
#| 
#| 
#| write_xlsx(
#|   mantel_FDR,
#|   "Mantel_Order_environment_FDR005.xlsx"
#| )
#| 
#| 
#| ############################################################
#| # 21. CORRELAZIONE TRA VARIABILI AMBIENTALI
#| ############################################################
#| #
#| # Usiamo Spearman perché:
#| #
#| # - n piccolo
#| # - non assumiamo normalità
#| # - alcune variabili sono fortemente correlate
#| #
#| ############################################################
#| 
#| corr_env <- linkET::correlate(
#|   env_scaled,
#|   method = "spearman"
#| )
#| 
#| 
#| ############################################################
#| # 22. CONTROLLO CORRELAZIONI
#| ############################################################
#| 
#| corr_env
#| 
#| 
#| ############################################################
#| # 23. GRAFICO MANTEL STILE DJP
#| ############################################################
#| 
#| p <- qcorrplot(
#|   corr_env,
#|   type = "lower",
#|   diag = FALSE
#| ) +
#|   
#|   ##########################################################
#| # MATRICE CORRELAZIONI
#| ##########################################################
#| 
#| geom_tile(
#|   aes(
#|     fill = r
#|   )
#| ) +
#|   
#|   ##########################################################
#| # ARCHI NON SIGNIFICATIVI
#| #
#| # Sono mostrati solamente per gli Order che presentano
#| # almeno una correlazione significativa.
#| ##########################################################
#| 
#| geom_couple(
#|   
#|   data = mantel_nonsig,
#|   
#|   aes(
#|     size = rd
#|   ),
#|   
#|   colour = "grey80",
#|   
#|   alpha = 0.6,
#|   
#|   curvature = nice_curvature()
#| ) +
#|   
#|   ##########################################################
#| # ARCHI SIGNIFICATIVI
#| ##########################################################
#| 
#| geom_couple(
#|   
#|   data = mantel_sig,
#|   
#|   aes(
#|     colour = pd,
#|     size = rd
#|   ),
#|   
#|   curvature = nice_curvature()
#| ) +
#|   
#|   ##########################################################
#| # COLORE MATRICE SPEARMAN
#| ##########################################################
#| 
#| scale_fill_gradientn(
#|   
#|   colours = brewer.pal(
#|     11,
#|     "RdBu"
#|   ),
#|   
#|   limits = c(
#|     -1,
#|     1
#|   )
#| ) +
#|   
#|   ##########################################################
#| # SPESSORE ARCHI = MANTEL r
#| ##########################################################
#| 
#| scale_size_manual(
#|   
#|   values = c(
#|     
#|     "< 0.2" = 0.4,
#|     
#|     "0.2 - 0.4" = 0.8,
#|     
#|     ">= 0.4" = 1.6
#|   ),
#|   
#|   drop = FALSE
#| ) +
#|   
#|   ##########################################################
#| # COLORE ARCHI = p-value
#| ##########################################################
#| 
#| scale_colour_manual(
#|   
#|   values = brewer.pal(
#|     3,
#|     "Set1"
#|   ),
#|   
#|   drop = TRUE
#| ) +
#|   
#|   ##########################################################
#| # LEGGENDE
#| ##########################################################
#| 
#| guides(
#|   
#|   size = guide_legend(
#|     
#|     title = "Mantel's r",
#|     
#|     override.aes = list(
#|       colour = "grey35"
#|     ),
#|     
#|     order = 2
#|   ),
#|   
#|   colour = guide_legend(
#|     
#|     title = "Mantel's p",
#|     
#|     override.aes = list(
#|       size = 3
#|     ),
#|     
#|     order = 1
#|   ),
#|   
#|   fill = guide_colorbar(
#|     
#|     title = "Spearman's r",
#|     
#|     order = 3
#|   )
#| ) +
#|   
#|   ##########################################################
#| # TEMA
#| ##########################################################
#| 
#| theme_minimal() +
#|   
#|   theme(
#|     
#|     axis.text.x = element_text(
#|       size = 9,
#|       angle = 45,
#|       hjust = 1
#|     ),
#|     
#|     axis.text.y = element_text(
#|       size = 9
#|     ),
#|     
#|     legend.title = element_text(
#|       size = 9
#|     ),
#|     
#|     legend.text = element_text(
#|       size = 8
#|     )
#|   )
#| 
#| 
#| ############################################################
#| # 24. MOSTRA FIGURA
#| ############################################################
#| 
#| p
#| 
#| 
#| ############################################################
#| # 25. ESPORTAZIONE PDF
#| ############################################################
#| 
#| ggsave(
#|   filename = "Mantel_Order_environment_ALL.pdf",
#|   plot = p,
#|   width = 11,
#|   height = 8,
#|   units = "in"
#| )
#| 
#| 
#| ############################################################
#| # 26. ESPORTAZIONE SVG
#| ############################################################
#| 
#| ggsave(
#|   filename = "Mantel_Order_environment_ALL.svg",
#|   plot = p,
#|   width = 11,
#|   height = 8,
#|   units = "in"
#| )
#| 
#| 
#| ############################################################
#| # 27. ESPORTAZIONE PNG ALTA RISOLUZIONE
#| ############################################################
#| 
#| ggsave(
#|   filename = "Mantel_Order_environment_ALL.png",
#|   plot = p,
#|   width = 11,
#|   height = 8,
#|   units = "in",
#|   dpi = 600
#| )
#| 
#| 
#| ############################################################
#| # 28. RIASSUNTO FINALE
#| ############################################################
#| 
#| cat(
#|   "\n============================================\n"
#| )
#| 
#| cat(
#|   "Numero campioni:",
#|   nrow(df_sel),
#|   "\n"
#| )
#| 
#| cat(
#|   "Numero Order testati:",
#|   ncol(df_sel),
#|   "\n"
#| )
#| 
#| cat(
#|   "Numero variabili ambientali:",
#|   ncol(env_scaled),
#|   "\n"
#| )
#| 
#| cat(
#|   "Order con almeno una associazione p < 0.05:",
#|   length(taxa_sig),
#|   "\n"
#| )
#| 
#| cat(
#|   "Associazioni totali p < 0.05:",
#|   nrow(mantel_sig),
#|   "\n"
#| )
#| 
#| cat(
#|   "Associazioni FDR < 0.05:",
#|   nrow(mantel_FDR),
#|   "\n"
#| )
#| 
#| cat(
#|   "============================================\n"
#| )
#| #################################################################################
#| mantel_significant_table %>%
#|   select(spec, env, r, p, FDR) %>%
#|   arrange(p)
#| length(taxa_sig)
#| taxa_sig
#| 
#| mantel <- mantel %>%
#|   ungroup() %>%
#|   mutate(
#|     FDR_global = p.adjust(p, method = "BH")
#|   )
#| 
#| mantel %>%
#|   filter(FDR < 0.05) %>%
#|   select(spec, env, r, p, FDR) %>%
#|   arrange(env, FDR, desc(r)) %>%
#|   print(n = Inf)
#| 
#| mantel %>%
#|   filter(FDR < 0.05) %>%
#|   distinct(spec) %>%
#|   pull(spec)
#| ################################################################################
#| ############################################################
#| # FIGURA MANTEL AGGIORNATA
#| # Solo associazioni significative dopo FDR < 0.05
#| ############################################################
#| 
#| library(tidyverse)
#| library(linkET)
#| library(RColorBrewer)
#| 
#| 
#| ############################################################
#| # 1. TAXA CON ALMENO UNA ASSOCIAZIONE FDR < 0.05
#| ############################################################
#| 
#| taxa_FDR <- mantel %>%
#|   filter(FDR < 0.05) %>%
#|   pull(spec) %>%
#|   unique()
#| 
#| 
#| cat(
#|   "\nTaxa con almeno una associazione FDR < 0.05:\n"
#| )
#| 
#| print(taxa_FDR)
#| 
#| cat(
#|   "\nNumero taxa selezionati:",
#|   length(taxa_FDR),
#|   "\n"
#| )
#| 
#| 
#| ############################################################
#| # 2. SOLO ASSOCIAZIONI SIGNIFICATIVE DOPO FDR
#| ############################################################
#| 
#| mantel_sig_FDR <- mantel %>%
#|   filter(
#|     spec %in% taxa_FDR,
#|     FDR < 0.05
#|   )
#| 
#| 
#| ############################################################
#| # 3. CONTROLLO RISULTATI
#| ############################################################
#| 
#| mantel_sig_FDR %>%
#|   select(
#|     spec,
#|     env,
#|     r,
#|     p,
#|     FDR,
#|     rd,
#|     pd
#|   ) %>%
#|   arrange(
#|     env,
#|     FDR,
#|     desc(r)
#|   ) %>%
#|   print(n = Inf)
#| 
#| 
#| ############################################################
#| # 4. SALVATAGGIO TABELLA DELLE ASSOCIAZIONI MOSTRATE
#| ############################################################
#| 
#| writexl::write_xlsx(
#|   mantel_sig_FDR,
#|   "Mantel_associations_FDR005_used_in_figure.xlsx"
#| )
#| 
#| 
#| ############################################################
#| # 5. CORRELAZIONI TRA VARIABILI AMBIENTALI
#| ############################################################
#| 
#| corr_env <- linkET::correlate(
#|   env_scaled,
#|   method = "spearman"
#| )
#| 
#| 
#| ############################################################
#| # 6. FIGURA
#| ############################################################
#| 
#| p_FDR <- qcorrplot(
#|   corr_env,
#|   type = "lower",
#|   diag = FALSE
#| ) +
#|   
#|   ##########################################################
#| # MATRICE DI CORRELAZIONE AMBIENTALE
#| ##########################################################
#| 
#| geom_tile(
#|   aes(
#|     fill = r
#|   )
#| ) +
#|   
#|   ##########################################################
#| # ARCHI MANTEL
#| # Solo FDR < 0.05
#| ##########################################################
#| 
#| geom_couple(
#|   data = mantel_sig_FDR,
#|   
#|   aes(
#|     colour = pd,
#|     size = rd
#|   ),
#|   
#|   curvature = nice_curvature()
#| ) +
#|   
#|   ##########################################################
#| # COLORI MATRICE SPEARMAN
#| ##########################################################
#| 
#| scale_fill_gradientn(
#|   
#|   colours = brewer.pal(
#|     11,
#|     "RdBu"
#|   ),
#|   
#|   limits = c(
#|     -1,
#|     1
#|   ),
#|   
#|   name = "Spearman's r"
#| ) +
#|   
#|   ##########################################################
#| # SPESSORE ARCHI = MANTEL r
#| ##########################################################
#| 
#| scale_size_manual(
#|   
#|   values = c(
#|     
#|     "< 0.2" = 0.4,
#|     
#|     "0.2 - 0.4" = 0.8,
#|     
#|     ">= 0.4" = 1.6
#|   ),
#|   
#|   drop = FALSE,
#|   
#|   name = "Mantel's r"
#| ) +
#|   
#|   ##########################################################
#| # COLORE ARCHI = p-value NOMINALE
#| ##########################################################
#| 
#| scale_colour_manual(
#|   
#|   values = c(
#|     "< 0.01" = brewer.pal(3, "Set1")[1],
#|     "0.01 - 0.05" = brewer.pal(3, "Set1")[2],
#|     ">= 0.05" = brewer.pal(3, "Set1")[3]
#|   ),
#|   
#|   drop = TRUE,
#|   
#|   name = "Mantel's p"
#| ) +
#|   
#|   ##########################################################
#| # LEGGENDE
#| ##########################################################
#| 
#| guides(
#|   
#|   colour = guide_legend(
#|     title = "Mantel's p",
#|     override.aes = list(
#|       size = 2
#|     ),
#|     order = 1
#|   ),
#|   
#|   size = guide_legend(
#|     title = "Mantel's r",
#|     override.aes = list(
#|       colour = "grey30"
#|     ),
#|     order = 2
#|   ),
#|   
#|   fill = guide_colorbar(
#|     title = "Spearman's r",
#|     order = 3
#|   )
#| ) +
#|   
#|   ##########################################################
#| # TEMA
#| ##########################################################
#| 
#| theme_minimal() +
#|   
#|   theme(
#|     
#|     axis.text.x = element_text(
#|       size = 9,
#|       angle = 45,
#|       hjust = 1
#|     ),
#|     
#|     axis.text.y = element_text(
#|       size = 9
#|     ),
#|     
#|     legend.title = element_text(
#|       size = 9
#|     ),
#|     
#|     legend.text = element_text(
#|       size = 8
#|     ),
#|     
#|     plot.margin = margin(
#|       10,
#|       25,
#|       10,
#|       10
#|     )
#|   )
#| 
#| 
#| ############################################################
#| # 7. MOSTRA FIGURA
#| ############################################################
#| 
#| p_FDR
#| 
#| 
#| ############################################################
#| # 8. SALVA PDF
#| ############################################################
#| 
#| ggsave(
#|   filename = "Mantel_Order_environment_FDR005.pdf",
#|   plot = p_FDR,
#|   width = 11,
#|   height = 8,
#|   units = "in"
#| )
#| 
#| 
#| ############################################################
#| # 9. SALVA SVG
#| ############################################################
#| 
#| ggsave(
#|   filename = "Mantel_Order_environment_FDR005.svg",
#|   plot = p_FDR,
#|   width = 11,
#|   height = 8,
#|   units = "in"
#| )
#| 
#| 
#| ############################################################
#| # 10. SALVA PNG 600 dpi
#| ############################################################
#| 
#| ggsave(
#|   filename = "Mantel_Order_environment_FDR005.png",
#|   plot = p_FDR,
#|   width = 11,
#|   height = 8,
#|   units = "in",
#|   dpi = 600
#| )
#| save.image()
#| 


