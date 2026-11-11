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


