# Load packages
library("dplyr", "DESeq2", "tximport", "GenomicFeatures")

# Set working directory containing subdirectories & sample names file
inDir <- "./"

# Read in sample names
samples <- read.table(file.path(inDir, "sample_names_ordered_DESeq2.txt"), header = FALSE)
colnames(samples) <- "names"
samples

# Define files
files <- file.path(inDir, samples$names, "quant.sf")
names(files) <- samples$names
all(file.exists(files))

# Import counts
countdata <- tximport(files, type = "salmon", txOut = TRUE)
names(countdata)
head(countdata$counts)

# Assign "treatment" and "genotype" variables.
(genotype <- factor(c(rep("WT", 18), rep("KO", 12))))
(treatment <- factor(c(rep("UT", 9), rep("IFNG_4h", 9), rep("UT", 6), rep("IFNG_4h", 6))))

# Create coldata table, matching sample name (row names) with treatment and replicate status
(coldata <- data.frame(row.names = colnames(countdata$counts), genotype, treatment))
coldata

# Instantiate the DESeqDataSet
dds <- DESeqDataSetFromTximport(countdata,
                                colData = coldata,
                                design = ~ genotype + treatment + genotype:treatment)

# Set the reference level
dds$treatment <- relevel(dds$treatment, ref = "UT")
dds$genotype <- relevel(dds$genotype, ref ="WT")

# Remove transcripts with zero counts
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds)

# Run DESeq2
dds <- DESeq(dds)
resultsNames(dds)

# A: Effect of treatment in WT (WT IFNG vs UT)
res_WT_IFNG_4h = results(dds, list(c("treatment_IFNG_4h_vs_UT")))
ix = which.min(res_WT_IFNG_4h$padj) # most significant
res_WT_IFNG_4h <- res_WT_IFNG_4h[order(res_WT_IFNG_4h$padj),] # sorts
table(res_WT_IFNG_4h$padj < 0.05)
head(res_WT_IFNG_4h)

# B: Effect of treatment in KOs (KO IFNG vs UT)
res_KO_IFNG_4h = results(dds, list(c("treatment_IFNG_4h_vs_UT", "genotypeKO.treatmentIFNG_4h")))
ix = which.min(res_KO_IFNG_4h$padj) # most significant
res_KO_IFNG_4h <- res_KO_IFNG_4h[order(res_KO_IFNG_4h$padj),] # sorts
table(res_KO_IFNG_4h$padj < 0.05)
head(res_KO_IFNG_4h)

# C: Difference between KO UT vs WT UT
res_KO_vs_WT_UT = results(dds, contrast = c("genotype", "KO", "WT"))
ix = which.min(res_KO_vs_WT_UT$padj) # most significant
res_KO_vs_WT_UT <- res_KO_vs_WT_UT[order(res_KO_vs_WT_UT$padj),] # sorts
table(res_KO_vs_WT_UT$padj < 0.05)
head(res_KO_vs_WT_UT)

# D: Difference between KO IFNG vs WT IFNG
res_KO_vs_WT_IFNG_4h = results(dds, list(c("genotype_KO_vs_WT","genotypeKO.treatmentIFNG_4h")))
ix = which.min(res_KO_vs_WT_IFNG_4h$padj) # most significant
res_KO_vs_WT_IFNG_4h <- res_KO_vs_WT_IFNG_4h[order(res_KO_vs_WT_IFNG_4h$padj),] # sorts
table(res_KO_vs_WT_IFNG_4h$padj < 0.05)
head(res_KO_vs_WT_IFNG_4h)

# E: Different responses in genotypes (interaction term) for KO
res_interaction = results(dds, name = "genotypeKO.treatmentIFNG_4h")
ix = which.min(res_interaction$padj) # most significant
res_interaction <- res_interaction[order(res_interaction$padj),] # sorts
table(res_interaction$padj < 0.05)
head(res_interaction)

# Extract normalized counts.
dds <- as.data.frame(counts(dds, normalized = TRUE))
colnames(dds)

# Merge results tables with normalized counts tables. Start with effect of treatment in WT
resdata_WT_IFNG_4h <- merge(as.data.frame(res_WT_IFNG_4h), as.data.frame(dds), by = "row.names", sort = FALSE)
names(resdata_WT_IFNG_4h)[1] <- "transcript"
head(resdata_WT_IFNG_4h)

# Effect of treatment in KO
resdata_KO_IFNG_4h <- merge(as.data.frame(res_KO_IFNG_4h), as.data.frame(dds), by = "row.names", sort = FALSE)
names(resdata_KO_IFNG_4h)[1] <- "transcript"
head(resdata_KO_IFNG_4h)

# Difference between KO and WT without treatment
resdata_KO_vs_WT_UT <- merge(as.data.frame(res_KO_vs_WT_UT), as.data.frame(dds), by = "row.names", sort = FALSE)
names(resdata_KO_vs_WT_UT)[1] <- "transcript"
head(resdata_KO_vs_WT_UT)

# Difference between KO and WT with treatment
resdata_KO_vs_WT_IFNG_4h <- merge(as.data.frame(res_KO_vs_WT_IFNG_4h), as.data.frame(dds), by = "row.names", sort = FALSE)
names(resdata_KO_vs_WT_IFNG_4h)[1] <- "transcript"
head(resdata_KO_vs_WT_IFNG_4h)

# Interaction term.
resdata_KO_interaction <- merge(as.data.frame(res_interaction), as.data.frame(dds), by = "row.names", sort = FALSE)
names(resdata_KO_interaction)[1] <- "transcript"
head(resdata_KO_interaction)

# Define outDir
outDir="./"

# Write output
write.table(resdata_WT_IFNG_4h, file.path(outDir, "DESeq2_WT_IFNG_4h_vs_UT.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(resdata_KO_IFNG_4h, file.path(outDir, "DESeq2_KO_IFNG_4h_vs_UT.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(resdata_KO_vs_WT_UT, file.path(outDir, "DESeq2_KO_vs_WT_UT.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(resdata_KO_vs_WT_IFNG_4h, file.path(outDir, "DESeq2_KO_vs_WT_IFNG_4h.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(resdata_KO_interaction, file.path(outDir, "DESeq2_interaction.txt"), sep = "\t", quote = FALSE, row.names = FALSE)