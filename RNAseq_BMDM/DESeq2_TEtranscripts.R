# Load packages
library("dplyr", "DESeq2", "apeglm")

# Set input directory
# This should be the directory that hosts your count table
workingDir <- "/Users/path/to/working/directory"

# Write session info
writeLines(capture.output(sessionInfo()), file.path(workingDir, "session_info.txt"))

# Read in counts table
countdata_picc_2h <- read.csv(file.path(workingDir, "piccolo_IFNG_2h_vs_UT_TEsonly.cntTable"), sep = "", header = FALSE)
countdata_plat_2h <- read.csv(file.path(workingDir, "platanitis_IFNG_2h_vs_UT_TEsonly.cntTable"), sep = "", header = FALSE)
countdata_picc_4h <- read.csv(file.path(workingDir, "piccolo_IFNG_4h_vs_UT_TEsonly.cntTable"), sep = "", header = FALSE)

# Rename columns
# TEtranscripts reorders the columns according to name
# Double check the ".cntTable" (not "TEsonly.cntTable") file for column order
colnames(countdata_picc_2h) <- c("family", "picc_IFNG_2h_R3", "picc_IFNG_2h_R1", "picc_IFNG_2h_R2", "picc_UT_R3", "picc_UT_R2", "picc_UT_R1")
colnames(countdata_picc_4h) <- c("family", "picc_IFNG_4h_R2", "picc_IFNG_4h_R1", "picc_IFNG_4h_R3", "picc_UT_R3", "picc_UT_R2", "picc_UT_R1")
colnames(countdata_picc_2h) <- c("family", "plat_IFNG_2h_R3", "plat_IFNG_2h_R2", "plat_IFNG_2h_R1", "plat_UT_R1", "plat_UT_R2", "plat_UT_R3")

# Optionally, reorder columns
countdata_picc_2h <- countdata_picc_2h[, c(1, 7, 6, 5, 3, 4, 2)]
countdata_picc_4h <- countdata_picc_4h[, c(1, 7, 6, 5, 3, 2, 4)]
countdata_plat_2h <- countdata_plat_2h[, c(1, 7, 6, 5, 4, 3, 2)]

# Set first column (family) to be the rownames
rownames(countdata_picc_2h) <- countdata_picc_2h[, which(colnames(countdata_picc_2h) == "family")]
countdata_picc_2h <- subset(countdata_picc_2h, select = -c(family))
rownames(countdata_picc_4h) <- countdata_picc_4h[, which(colnames(countdata_picc_4h) == "family")]
countdata_picc_4h <- subset(countdata_picc_4h, select = -c(family))
rownames(countdata_plat_2h) <- countdata_plat_2h[, which(colnames(countdata_plat_2h) == "family")]
countdata_plat_2h <- subset(countdata_plat_2h, select = -c(family))

# Filter out outliers (clustered separately from other samples by PCA pre- and post-DESeq2 analysis)
countdata_picc_2h <- select(countdata_picc_2h, -UT_R2)
countdata_picc_4h <- select(countdata_picc_4h, -UT_R2)
countdata_plat_2h <- select(countdata_plat_2h, -UT_R3)

# Only retain rows with at least one read.
countdata_picc_2h <- countdata_picc_2h[apply(countdata_picc_2h,1,function(x){max(x)}) > 1,]
countdata_picc_4h <- countdata_picc_4h[apply(countdata_picc_4h,1,function(x){max(x)}) > 1,]
countdata_plat_2h <- countdata_plat_2h[apply(countdata_plat_2h,1,function(x){max(x)}) > 1,]

# Convert countdata tables into a matrix - necessary for running DESeq2.
countdata_picc_2h <- as.matrix(countdata_picc_2h)
countdata_picc_4h <- as.matrix(countdata_picc_4h)
countdata_plat_2h <- as.matrix(countdata_plat_2h)

# Assign control vs treat samples
(treatment_picc_2h <- factor(c(rep("UT", 2), rep("IFNG_2h", 3))))
(treatment_picc_4h <- factor(c(rep("UT", 2), rep("IFNG_4h", 3))))
(treatment_plat_2h <- factor(c(rep("UT", 2), rep("IFNG_2h", 3))))

# Create a "coldata" table containing the sample names with their appropriate condition
(coldata_picc_2h <- data.frame(row.names = colnames(countdata_picc_2h), treatment_picc_2h))
(coldata_picc_4h <- data.frame(row.names = colnames(countdata_picc_4h), treatment_picc_4h))
(coldata_plat_2h <- data.frame(row.names = colnames(countdata_plat_2h), treatment_plat_2h))

# Construct a DESeqDataSet
dds_picc_2h <- DESeqDataSetFromMatrix(countData = countdata_picc_2h, colData = coldata_picc_2h, design = ~ treatment_picc_2h)
dds_picc_4h <- DESeqDataSetFromMatrix(countData = countdata_picc_4h, colData = coldata_picc_4h, design = ~ treatment_picc_4h)
dds_plat_2h <- DESeqDataSetFromMatrix(countData = countdata_plat_2h, colData = coldata_plat_2h, design = ~ treatment_plat_2h)

# Set the reference level
dds_picc_2h$treatment_picc_2h <- relevel(dds_picc_2h$treatment_picc_2h, ref = "UT")
dds_picc_4h$treatment_picc_4h <- relevel(dds_picc_4h$treatment_picc_4h, ref = "UT")
dds_plat_2h$treatment_plat_2h <- relevel(dds_plat_2h$treatment_plat_2h, ref = "UT")

# Run DESeq2
dds_picc_2h <- DESeq(dds_picc_2h)
dds_plat_2h <- DESeq(dds_plat_2h)
dds_plat_2h <- DESeq(dds_plat_2h)

# Get differential expression results
res_picc_2h <- results(dds_picc, contrast=c("treatment_picc", "IFNG_2h", "UT"))
res_picc_4h <- results(dds_picc, contrast=c("treatment_picc", "IFNG_4h", "UT"))
res_plat_2h <- results(dds_picc, contrast=c("treatment_plat", "IFNG_2h", "UT"))

# Get differential expression results, this time shrinking by apeglm for visualization
resLFC_picc_2h <- lfcShrink(dds_picc, coef="treatment_picc_IFNG_2h_vs_UT", type="apeglm")
resLFC_picc_4h <- lfcShrink(dds_picc, coef="treatment_picc_IFNG_4h_vs_UT", type="apeglm")
resLFC_plat_2h <- lfcShrink(dds_picc, coef="treatment_plat_IFNG_2h_vs_UT", type="apeglm")

# Extract normalized counts
dds_picc_2h <- as.data.frame(counts(dds_picc_2h, normalized = TRUE))
dds_picc_4h <- as.data.frame(counts(dds_picc_4h, normalized = TRUE))
dds_plat_2h <- as.data.frame(counts(dds_plat_2h, normalized = TRUE))

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

# Write DESeq2 output files
write.table(resdata_picc_2h, file.path(workingDir, "picc_BMDM_IFNG_2h_vs_UT_TEtranscripts.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_picc_4h, file.path(workingDir, "picc_BMDM_IFNG_4h_vs_UT_TEtranscripts.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_plat_2h, file.path(workingDir, "plat_BMDM_IFNG_2h_vs_UT_TEtranscripts.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdataLFC_picc_2h, file.path(workingDir, "picc_BMDM_IFNG_2h_vs_UT_LFC_apeglm_TEtranscripts.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdataLFC_picc_4h, file.path(workingDir, "picc_BMDM_IFNG_4h_vs_UT_LFC_apeglm_TEtranscripts.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdataLFC_plat_2h, file.path(workingDir, "plat_BMDM_IFNG_2h_vs_UT_LFC_apeglm_TEtranscripts.txt"), quote = FALSE, row.names = FALSE, sep = "\t")