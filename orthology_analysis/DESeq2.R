# Load packages
library("dplyr", "DESeq2", "apeglm")

# Set input directory
# This should be the directory that hosts your count table
workingDir <- "/Users/path/to/working/directory"

# Write session info
writeLines(capture.output(sessionInfo()), file.path(workingDir, "session_info.txt"))

# Read in full counts table
countdata <- read.csv(file.path(workingDir, "raw_gene_counts.txt"), sep = "", header = TRUE, skip = 1)
colnames(countdata)[colnames(countdata) == "Geneid"] <- "gene_id"

# Remove the first six columns (geneid, chr, start, end, strand, length)
countdata <- countdata[ ,-c(1:6)]

# Remove suffix and path (if necessary) from column names
colnames(countdata) <- gsub("\\.sorted.uniq.bam$", "", colnames(countdata))
colnames(countdata) <- gsub("\\X.Users.path.to.working.directory", "", colnames(countdata))
colnames(countdata) <- gsub("\\X.Users.path.to.working.directory", "", colnames(countdata))

# Set first column (gene_id) to be the rownames
rownames(countdata) <- countdata[, which(colnames(countdata) == "gene_id")]
countdata <- subset(countdata, select = -gene_id)

# Assign control vs treat samples
(treatment <- factor(c(rep("UT", 2), rep("IFNG_24h", 2))))

# Create a "coldata" table containing the sample names with their appropriate condition
(coldata <- data.frame(row.names = colnames(countdata), treatment))

# Construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ treatment)

# Set the reference level
dds$treatment <- relevel(dds$treatment, ref = "UT")

# Remove genes with zero counts
dds <- dds[rowSums(counts(dds)) > 1,]

# Run DESeq2
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds, contrast=c("treatment", "IFNG_24h", "UT"))

# Extract normalized counts
dds <- as.data.frame(counts(dds, normalized = TRUE))

# Merge unshrinked results table with normalized counts tables
resdata <- merge(as.data.frame(res), as.data.frame(dds), by = "row.names", sort = FALSE)
names(resdata)[1] <- "gene"

# Filter
resdata_padj <- subset(resdata_IFNG_24h, padj <= 0.05)
resdata_padj_log2FC0_up <- subset(resdata_padj, log2FoldChange >= 0)

# Write DESeq2 output file
write.table(resdata, file.path(workingDir, "qiao_mono_IFNG_24h.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_padj_log2FC0_up, file.path(workingDir, "qiao_mono_IFNG_24h_padj_log2FC0_up.txt"), quote = FALSE, row.names = FALSE, sep = "\t")