# Load required packages
# Credit for function: https://gist.github.com/smithdanielle/9913897
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}
packages <- c("dplyr", "DESeq2", "apeglm")
check.packages(packages)

# Set input directory
# This should be the directory that hosts your count table
workingDir <- "/Users/path/to/working/directory"

# Write session info
writeLines(capture.output(sessionInfo()), file.path(workingDir, "session_info.txt"))

# Read in full counts table
countdata <- read.csv(file.path(workingDir, "raw_gene_counts.txt"), sep = "", header = TRUE, skip = 1)

# Read in crossref table
# This provides the gene name for each gene ID
crossref <- read.csv(file.path(workingDir, "gencode_vM18_crossref.txt"), sep = "", header = FALSE)
colnames(crossref) <- c("gene_id", "gene")

# Remove the first six columns (geneid, chr, start, end, strand, length)
countdata <- countdata[ ,-c(1:6)]

# Remove suffix and path (if necessary) from column names
colnames(countdata) <- gsub("\\.sorted.uniq.bam$", "", colnames(countdata))
colnames(countdata) <- gsub("\\X.Users.path.to.working.directory", "", colnames(countdata))
colnames(countdata) <- gsub("\\X.Users.path.to.working.directory", "", colnames(countdata))

# Match gene names to gene IDs
countdata <- left_join(countdata, crossref, by = "gene_id")

# Set first column (geneid) to be the rownames
rownames(countdata) <- countdata[, which(colnames(countdata) == "gene")]
countdata <- subset(countdata, select = -c(gene, gene_name))

# Parse table by dataset
countdata_picc <- select(countdata, contains("picc_"))
countdata_plat <- select(countdata, contains("plat_"))

# Convert countdata tables into a matrix
countdata_picc <- as.matrix(countdata_picc)
countdata_plat <- as.matrix(countdata_plat)

# Assign control vs treat samples
(treatment_picc <- factor(c(rep("UT", 3), rep("IFNG_2h", 3), rep("IFNG_4h", 3))))
(treatment_plat <- factor(c(rep("UT", 3), rep("IFNG_2h", 3))))

# Create a "coldata" table containing the sample names with their appropriate condition
(coldata_picc <- data.frame(row.names = colnames(countdata_picc), treatment_picc))
(coldata_plat <- data.frame(row.names = colnames(countdata_plat), treatment_plat))

# Construct a DESeqDataSet
dds_picc <- DESeqDataSetFromMatrix(countData = countdata_picc, colData = coldata_picc, design = ~ treatment_picc)
dds_plat <- DESeqDataSetFromMatrix(countData = countdata_plat, colData = coldata_plat, design = ~ treatment_plat)

# Set the reference level
dds_picc$treatment_picc <- relevel(dds_picc$treatment_picc, ref = "UT")
dds_plat$treatment_plat <- relevel(dds_plat$treatment_plat, ref = "UT")

# Remove genes with zero counts
dds_picc <- dds_picc[rowSums(counts(dds_picc)) > 1,]
dds_plat <- dds_plat[rowSums(counts(dds_plat)) > 1,]

# Run DESeq2
dds_picc <- DESeq(dds_picc)
dds_plat <- DESeq(dds_plat)

# Get differential expression results
res_picc_2h <- results(dds_picc, contrast=c("treatment_picc", "IFNG_2h", "UT"))
res_picc_4h <- results(dds_picc, contrast=c("treatment_picc", "IFNG_4h", "UT"))
res_plat_2h <- results(dds_picc, contrast=c("treatment_plat", "IFNG_2h", "UT"))

# Get differential expression results, this time shrinking by apeglm for visualization
resLFC_picc_2h <- lfcShrink(dds_picc, coef="treatment_picc_IFNG_2h_vs_UT", type="apeglm")
resLFC_picc_4h <- lfcShrink(dds_picc, coef="treatment_picc_IFNG_4h_vs_UT", type="apeglm")
resLFC_plat_2h <- lfcShrink(dds_picc, coef="treatment_plat_IFNG_2h_vs_UT", type="apeglm")

# Extract normalized counts
dds_picc <- as.data.frame(counts(dds_picc, normalized = TRUE))
dds_plat <- as.data.frame(counts(dds_plat, normalized = TRUE))

# Merge unshrinked results table with normalized counts tables
resdata_picc_2h <- merge(as.data.frame(res_picc_2h), as.data.frame(dds_picc), by = "row.names", sort = FALSE)
names(resdata_picc_2h)[1] <- "gene"
resdata_picc_4h <- merge(as.data.frame(res_picc_4h), as.data.frame(dds_picc), by = "row.names", sort = FALSE)
names(resdata_picc_4h)[1] <- "gene"
resdata_plat_2h <- merge(as.data.frame(res_plat_2h), as.data.frame(dds_plat), by = "row.names", sort = FALSE)
names(resdata_plat_2h)[1] <- "gene"

# Merge shrinked results table with normalized counts tables
resdataLFC_picc_2h <- merge(as.data.frame(resLFC_picc_2h), as.data.frame(dds_picc), by = "row.names", sort = FALSE)
names(resdataLFC_picc_2h)[1] <- "gene"
resdataLFC_picc_4h <- merge(as.data.frame(resLFC_picc_4h), as.data.frame(dds_picc), by = "row.names", sort = FALSE)
names(resdataLFC_picc_4h)[1] <- "gene"
resdataLFC_plat_2h <- merge(as.data.frame(resLFC_plat_2h), as.data.frame(dds_plat), by = "row.names", sort = FALSE)
names(resdataLFC_plat_2h)[1] <- "gene"

# Subset for genes w/ padj >= 0.90, baseMean >= 100, and abs(log2FC) <= 0.10
resdata_nonresponsive_picc_2h <- filter(resdata_picc_2h, baseMean >= 100, padj >= 0.90, abs(log2FoldChange) <= 0.10)
resdata_nonresponsive_picc_4h <- filter(resdata_picc_4h, baseMean >= 100, padj >= 0.90, abs(log2FoldChange) <= 0.10)
resdata_nonresponsive_plat_2h <- filter(resdata_plat_2h, baseMean >= 100, padj >= 0.90, abs(log2FoldChange) <= 0.10)

# Read in bed file that provides TSS coordinates for each gene
all_gene_tss <- read.table(file.path(workingDir, "gencode_vM18_tss.bed"), sep="\t", header = FALSE)
colnames(all_gene_tss) <- c("chr", "start", "end", "gene_id", "gene", "strand")
all_gene_tss <- select(all_gene_tss, c(1:4, 6))

# Append unique gene names using crossref file
all_gene_tss <- right_join(all_gene_tss, crossref, by = "gene_id") %>% select(c(1:4,6,5))

# Make bed files
bed_picc_2h <- right_join(all_gene_tss, resdata_picc_2h, by = "gene") %>% select(1:6)
bed_picc_4h <- right_join(all_gene_tss, resdata_picc_4h, by = "gene") %>% select(1:6)
bed_plat_2h <- right_join(all_gene_tss, resdata_plat_2h, by = "gene") %>% select(1:6)
bed_nonresponsive_picc_2h <- right_join(all_gene_tss, resdata_nonresponsive_picc_2h, by = "gene") %>% select(1:6)
bed_nonresponsive_picc_4h <- right_join(all_gene_tss, resdata_nonresponsive_picc_4h, by = "gene") %>% select(1:6)
bed_nonresponsive_plat_2h <- right_join(all_gene_tss, resdata_nonresponsive_plat_2h, by = "gene") %>% select(1:6)

# Write DESeq2 output files
write.table(resdata_picc_2h, file.path(workingDir, "picc_BMDM_IFNG_2h_vs_UT.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_picc_4h, file.path(workingDir, "picc_BMDM_IFNG_4h_vs_UT.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_plat_2h, file.path(workingDir, "plat_BMDM_IFNG_2h_vs_UT.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdataLFC_picc_2h, file.path(workingDir, "picc_BMDM_IFNG_vs_UT_2h_LFC_apeglm.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdataLFC_picc_4h, file.path(workingDir, "picc_BMDM_IFNG_vs_UT_4h_LFC_apeglm.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdataLFC_plat_2h, file.path(workingDir, "plat_BMDM_IFNG_vs_UT_2h_LFC_apeglm.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_nonresponsive_picc_2h, file.path(workingDir, "picc_BMDM_IFNG_2h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10.txt"),  quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_nonresponsive_picc_4h, file.path(workingDir, "picc_BMDM_IFNG_4h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10.txt"),  quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_nonresponsive_plat_2h, file.path(workingDir, "plat_BMDM_IFNG_2h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10.txt"),  quote = FALSE, row.names = FALSE, sep = "\t")

# Write bed output files
write.table(bed_picc_2h, file.path(workingDir, "picc_BMDM_IFNG_2h_vs_UT.bed"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(bed_picc_4h, file.path(workingDir, "picc_BMDM_IFNG_4h_vs_UT.bed"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(bed_plat_2h, file.path(workingDir, "plat_BMDM_IFNG_2h_vs_UT.bed"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(bed_nonresponsive_picc_2h, file.path(workingDir, "picc_BMDM_IFNG_2h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10.bed"),  quote = FALSE, row.names = FALSE, sep = "\t")
write.table(bed_nonresponsive_picc_4h, file.path(workingDir, "picc_BMDM_IFNG_4h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10.bed"),  quote = FALSE, row.names = FALSE, sep = "\t")
write.table(bed_nonresponsive_plat_2h, file.path(workingDir, "plat_BMDM_IFNG_2h_vs_UT_nonresponsive_padj0.90_baseMean100_log2FC0.10.bed"),  quote = FALSE, row.names = FALSE, sep = "\t")