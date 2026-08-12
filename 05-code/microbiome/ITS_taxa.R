# .libPaths('/home/szprisco/shared/MyRlibs_ngo00103') # this is the path file of where MSI locally installed dada2- might change if they have to reinstall
# setwd("/scratch.global/wupchurc/ITS2/02_filtered") # this is an example, input the path file of where the files are located

library(dada2)
library(phyloseq)

path <- ("/scratch.global/wupchurc/Blumeria-ITS2/02_filtered") # put the directory path of where the ITS2 sequencing files are-- they have L001 in the file name!
list.files(path) # check that all of the ITS2 files are there

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) # forward reads
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE)) # reverse reads
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "R1"), `[`, 1) # this is important to keep track of sample names

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# filter and trim raw reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     maxN = 0, maxEE = c(2, 4), truncQ = 2, 
                     minLen = 50, rm.phix = TRUE, 
                     compress = TRUE, multithread = FALSE)
head(out)

#dereplicate reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

# error models
errF <- learnErrors(derep_forward, multithread=8, randomize=TRUE)
errR <- learnErrors(derep_reverse, multithread=8, randomize=TRUE)

dadaFs <- dada(derep_forward, err=errF, multithread=8, pool="pseudo")
dadaRs <- dada(derep_reverse, err=errR, multithread=8, pool="pseudo")

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, 
                               derep_reverse, trimOverhang=TRUE, minOverlap=20)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, "/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/03-analysis/seqtab.rds") # this saves the output into a variable file, so if you can't finish the rest of the analysis the output file is saved into a file. When you come back, you can call this RDS file and continue the analysis

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=8, verbose=TRUE)
dim(seqtab.nochim)
lname <- nchar(colnames(seqtab.nochim))
summary(lname)
# seqtab.nochim <- seqtab.nochim[,(lname > 280)]
saveRDS(seqtab.nochim, "/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/03-analysis/seqtab_nochim.rds") # this no chimera file is the one that will be used, so make sure to save this rds file if you run out of time

#seqtab.nochim <- readRDS("/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/seqtab_nochim.rds")

uniquesToFasta(seqtab.nochim, fout = "/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/03-analysis/sequences.fasta")

# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                          filtered=out[,2], dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN),nonchim=rowSums(seqtab.nochim))
write.table(summary_tab, file = "/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/03-analysis/sequence_process_summary.txt", sep = "\t", quote=FALSE)

ref <- "/common/bioref/microbiome/dada2_taxonomy_references/sh_general_release_dynamic_all_19.02.2025.fasta"
# ref <- "/scratch.global/wupchurc/unite_2025/sh_general_release_dynamic_19.02.2025.fasta"
# ref <- "/common/bioref/microbiome/dada2_taxonomy_references/sh_general_release_dynamic_all_04.02.2020.fasta"
# ref <- "/common/bioref/microbiome/dada2_taxonomy_references/sh_general_release_dynamic_s_all_25.07.2023_dev.fasta"
# ref <- "/common/bioref/microbiome/dada2_taxonomy_references/sh_general_release_dynamic_s_all_04.02.2020.fasta"

# taxa <- assignTaxonomy(seqtab.nochim, ref, multithread = TRUE, outputBootstraps = TRUE, tryRC = TRUE)
taxa <- assignTaxonomy(seqtab.nochim, ref, multithread = TRUE, outputBootstraps = TRUE)

taxout <- taxa$tax
bootout <- taxa$boot


# both1 <- cbind(t(seqtab.nochim),taxa$tax, taxa$boot)
# write.table(both1, 
            # file = "/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/ITS_combined_sequences_taxa_bootstrap_2020ref_nofilter.txt",
            # sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxa$tax)
write.table(both2, 
            file = "/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/04-results/ITS2_taxa.txt", 
            sep = "\t", quote = FALSE, col.names=NA)


# ==============================================================================
# PHYLOSEQ: GENUS-LEVEL AGGLOMERATION & RELATIVE ABUNDANCE
# ==============================================================================

# 1. Build the phyloseq object
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxout)
)

# 2. Agglomerate data at the Genus level
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

# 3. Transform to Relative Abundance (proportions summing to 1 per sample)
# Note: If you want percentages (0-100%), change 'x / sum(x)' to '100 * x / sum(x)'
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# 4. Export the total Genus-level relative abundance table
genus_counts <- t(as(otu_table(ps_genus_rel), "matrix"))
genus_taxonomy <- as(tax_table(ps_genus_rel), "matrix")

genus_export <- cbind(genus_counts, genus_taxonomy)

write.table(genus_export, 
            file = "/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/04-results/ITS2_taxa_genus_relative_abundance.txt", 
            sep = "\t", quote = FALSE, col.names = NA)


# ==============================================================================
# FUNGI ONLY EXPORT (RELATIVE ABUNDANCE)
# ==============================================================================

# 1. Filter to Fungi only
ps_f <- subset_taxa(ps, Kingdom == "k__Fungi")

# 2. Agglomerate Fungi at the Genus level
ps_f_genus <- tax_glom(ps_f, taxrank = "Genus", NArm = FALSE)

# 3. Transform Fungi to Relative Abundance
# SEE NOTE BELOW regarding the ordering of subsetting vs transforming!
ps_f_genus_rel <- transform_sample_counts(ps_f_genus, function(x) x / sum(x))

# 4. Export the Fungi Genus-level relative abundance table
genus_counts_f <- t(as(otu_table(ps_f_genus_rel), "matrix"))
genus_taxonomy_f <- as(tax_table(ps_f_genus_rel), "matrix")

genus_export_f <- cbind(genus_counts_f, genus_taxonomy_f)

write.table(genus_export_f, 
            file = "/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/04-results/ITS_taxa_genus_relative_abundance_fungi_only.txt", 
            sep = "\t", quote = FALSE, col.names = NA)

