---
title: "deseq2"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Read in featureCounts output and run DESeq2
```{r}
## Set the directroy containing counts table
dir <- "/Users/%u/path_to_featureCounts_table"
## Read in counts table
countdata_IDs <- read.csv(file.path(dir, "featureCounts.txt"), sep = "", header = TRUE, skip = 1)
head(countdata_IDs)
```

Read in crossref table
```{r}
## Set the directory containing crossref table
dir2 <- "/Shares/CL_Shared/db/genomes/mm10/annotations"

## Read in crossref table
crossref <- read.csv(file.path(dir2, "geneSymbol_geneName_noBlanks.txt"), sep = "", header = FALSE)
colnames(crossref) <- c("gene", "gene_name") #changes column names to "gene" (col 1 has gene IDs) and "gene name" (col 2 has gene names)
head(crossref)
```

Filter out unnecessary columns (chr, start, end, strand, and length) & edit column names
```{r}
## Filter out unnecessary columns (chr, start, end, strand, and length) & rename "Geneid" to "gene"
countdata_IDs <- subset(countdata_IDs, select = -c(Chr, Start, End, Strand, Length)) #removes all of these columns
colnames(countdata_IDs)[colnames(countdata_IDs) == "Geneid"] <- "gene" #renames Geneid column to "gene"
head(countdata_IDs)

## Edit the column names
colnames(countdata_IDs) <- gsub("\\.uniq.bam$", "", colnames(countdata_IDs))
colnames(countdata_IDs) <- gsub("\\X.Users.%u.", "",
                               colnames(countdata_IDs))
colnames(countdata_IDs)
head(countdata_IDs)
```

```{r}
## Load dplyr
library(dplyr)

## Match gene names to gene IDs
countdata <- left_join(countdata_IDs, crossref, by = "gene") #adds gene_name column to end of table, creates new table now labeled "countdata_names", joins rows from crossref that match rows in countdata into countdata table
head(countdata)

## Discard gene ID column
countdata <- select(countdata, -gene) #deletes column labeled "gene" which contains IDs
head(countdata)

## Set gene names as row names & remove gene name column
rownames(countdata) <- countdata$gene_name #put gene names into rows
countdata <- select(countdata, -gene_name) #deletes column with gene names
head(countdata)
```

```{r}
## Assign "treatment" variables - i.e. duration of stimulation
(treatment <- factor(c("KO1","KO1","KO1","WT","WT","WT")))
(replicate <- factor(c(1,2,3,1,2,3))) 
```

```{r}
library("DESeq2")
## Create coldata table, matching sample name (row names) with treatment and replicate status
(coldata <- data.frame(row.names = colnames(countdata), treatment, replicate))
## Instantiate the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ treatment)
## Set the reference level
dds$treatment <- relevel(dds$treatment, ref = "WT")
```

## Run the DESeq pipeline.
```{r}
## Remove genes with zero counts
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds)
## Run DESeq2
dds <- DESeq(dds)
resultsNames(dds)
```

## Plot dispersion estimates.
```{r}
## Set output directory for writing png files
pngDir <- "/Users/%u/path_to_png"
# Plot dispersion estimates
plotDispEsts(dds)
dev.copy(png, res=200, height = 1500, width = 1500, pointsize=10, file.path(pngDir, "dispEst.png"))
dev.off()
```

## Get differential expression results and list all entries that have an adjusted p-value less than 0.05.
```{r}
res_KO1 <- results(dds, contrast=c("treatment", "KO1", "WT"))
res_KO1 <- res_KO1[order(res_KO1$padj),]
table(res_KO1$padj < 0.05)

#Log2FC shrinkage table
library("apeglm")
resLFC_KO1 <- lfcShrink(dds, coef="treatment_KO1_vs_WT", type="apeglm")
resLFC_KO1 <- resLFC_KO1[order(resLFC_KO1$padj),] #sorts table by padj value
table(resLFC_KO1$padj < 0.05)
```

## Make histograms of adjusted p values.
```{r}
hist(res_KO1$padj[res_KO1$baseMean > 1], 
     breaks = 0:20/20, 
     col = "grey50", 
     border = "white", 
     ylim = c(0, 2000), 
     xlab = "padj for genes with mean normalized count > 1", 
     main = "Histogram of padj for genes with mean normalized count > 1")
dev.copy(png, res=200, height = 1500, width = 1500, pointsize=10, file.path(pngDir, "hist_padj.png"))
dev.off()
```

## Make MA plots.
```{r}
plotMA(res_KO1, ylim = c(-10, 10), alpha = 0.05)
dev.copy(png, res = 200, height = 1500, width = 1500, pointsize = 10, file.path(pngDir, "MA.png"))
dev.off()
```

## Calculate sample distance from count data transformed using VST & visualize as a heatmap
```{r}
library("vsn")
library("dplyr")
library("ggplot2")
library("hexbin")
library("pheatmap")
library("RColorBrewer")
## Transform count data via VST
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)
meanSdPlot(assay(vsd), ranks = FALSE)
## Generate scatterplot using VST
df <- bind_rows(as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst")) 
colnames(df)[1:2] <- c("x", "y")
ggplot(df, aes(x = x, y = y)) + 
  geom_hex(bins = 80) +
  coord_fixed() + 
  facet_grid(. ~ transformation)
## Calculate distance
sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( "J774", vsd$treatment, vsd$replicate, sep = "_") 
colnames(sampleDistMatrix) <- paste( "J774", vsd$treatment, vsd$replicate, sep = "_")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.copy(png, res=200, height = 1000, width = 2000, pointsize=4, file.path(pngDir, "heatmap_VST.png"))
dev.off()
```

## Calculate sample distance from count data transformed using rlog & visualize as a heatmap
```{r}
# Transform count data via rlog 
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
meanSdPlot(assay(rld), ranks = FALSE)
## Generate scatterplot using rlog
df <- bind_rows(
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df)[1:2] <- c("x", "y")
ggplot(df, aes(x = x, y = y)) + 
  geom_hex(bins = 80) +
  coord_fixed() + 
  facet_grid( . ~ transformation)
## Calculate distance
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists)
rownames(sampleDistMatrix) <- paste( "J774", rld$treatment, rld$replicate, sep = "_") 
colnames(sampleDistMatrix) <- paste( "J774", rld$treatment, rld$replicate, sep = "_")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.copy(png, res=200, height = 1000, width = 2000, pointsize=4, file.path(pngDir, "heatmap_rlog.png"))
dev.off()
```

## Calculate sample distances using Poisson distance & make a heatmap
```{r}
library("PoiClaClu")
library("pheatmap")
library("RColorBrewer")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd)
rownames(samplePoisDistMatrix) <- paste( "J774", dds$treatment, dds$replicate, sep="_") 
colnames(samplePoisDistMatrix) <- paste( "J774", dds$treatment, dds$replicate, sep="_")
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
dev.copy(png, res=200, height = 1000, width = 2000, pointsize=4, file.path(pngDir, "heatmap_poisson.png"))
dev.off()
```

## Visualize VST distances in a PCA plot.
```{r}
library("affy")
pcaData <- plotPCA(vsd, intgroup = c("treatment"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"), digits = 2)
ggplot(pcaData, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size =3) +
  ylim(-2.5,3.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 3) + 
  labs(title = "J774 IFNG-treated WT (pooled) vs KO1", subtitle = "PCA - VST transformed")
dev.copy(png, res=200, height = 1000, width = 1000, pointsize=4, file.path(pngDir, "PCA_VST.png"))
dev.off()
```

## Visualize rlog distances in a PCA plot.
```{r}
library("affy")
pcaData <- plotPCA(rld, intgroup = c( "treatment"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"), digits = 2)
ggplot(pcaData, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size =3) +
  ylim(-3.6,5.6) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 2) + 
  labs(title = "J774 IFNG-treated WT (pooled) vs KO1", subtitle = "PCA - rlog transformed")
dev.copy(png, res=200, height = 1000, width = 1000, pointsize=4, file.path(pngDir, "PCA_rlog.png"))
dev.off()
```

## Generate heatmaps showing hierarchal clustering of top 50 or 20 genes with the highest variance.
```{r}
library("genefilter")
library("data.table")
library("dplyr")
## Get top 100 genes showing greatest variability & make heatmap
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100)
mat  <- assay(vsd)[topVarGenes,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, "treatment"])
rownames(anno) <- colnames(mat)
colnames(anno)[1] <- "treatment group"
pheatmap(mat, annotation_col = anno, fontsize_row = 12, fontsize_col = 12)
dev.copy(png, res=200, height = 6000, width = 3000, pointsize=4, file.path(pngDir, "geneClusterHeatmap_VST_topVar100.png"))
dev.off()
## Get top 50 genes showing greatest variability & make heatmap
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat  <- assay(vsd)[topVarGenes,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, "treatment"])
rownames(anno) <- colnames(mat)
colnames(anno)[1] <- "treatment group"
pheatmap(mat, annotation_col = anno, fontsize_row = 12, fontsize_col = 12)
dev.copy(png, res=200, height = 3000, width = 1500, pointsize=4, file.path(pngDir, "geneClusterHeatmap_VST_topVar50.png"))
dev.off()
## Get top 20 genes showing greatest variability & make heatmap
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[topVarGenes,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, "treatment"])
rownames(anno) <- colnames(mat)
colnames(anno)[1] <- "treatment group"
pheatmap(mat, annotation_col = anno, fontsize_row = 12, fontsize_col = 12)
dev.copy(png, res=200, height = 1300, width = 1300, pointsize=4, file.path(pngDir, "geneClusterHeatmap_VST_topVar20.png"))
dev.off()
```

## Merge the results tables with the normalized count data.
```{r}
## Extract normalized counts
dds_KO1 <- as.data.frame(counts(dds, normalized = TRUE))
colnames(dds_KO1)
## Merge results tables with normalized counts tables
resdata_KO1 <- merge(as.data.frame(res_KO1), as.data.frame(dds_KO1), by = "row.names", sort = FALSE)
names(resdata_KO1)[1] <- "gene"
head(resdata_KO1)
```

## Filter out entires that have an adjusted p-value greater than 0.05.
```{r}
nrow(resdata_KO1)
resdata_KO1_padj <- subset(resdata_KO1, padj <= 0.05)
nrow(resdata_KO1_padj)
```

## Subset by FC such that you have two tables per sample, one containing genes FC > 1 and another containing genes FC < -1.
```{r}
## Subset for genes w/ log2FC > 0 (less stringent)
nrow(resdata_KO1_padj)
resdata_KO1_padj_up_log2FC0 <- subset(resdata_KO1_padj, log2FoldChange >= 0)
nrow(resdata_KO1_padj_up_log2FC0)
## Subset for genes w/ log2FC > 1 (more stringent)
nrow(resdata_KO1_padj)
resdata_KO1_padj_up_log2FC1 <- subset(resdata_KO1_padj, log2FoldChange >= 1)
nrow(resdata_KO1_padj_up_log2FC1)
## Subset for genes that are downregulated (log2FC less than -1)
nrow(resdata_KO1_padj)
resdata_KO1_padj_down_log2FC0 <- subset(resdata_KO1_padj, log2FoldChange <= 0)
nrow(resdata_KO1_padj_down_log2FC0)
## Subset for genes that are downregulated (log2FC less than -1)
nrow(resdata_KO1_padj)
resdata_KO1_padj_down_log2FC1 <- subset(resdata_KO1_padj, log2FoldChange <= -1)
nrow(resdata_KO1_padj_down_log2FC1)
```

## Write results.
```{r}
## Set output directory
outDir <- "/Users/%u/path_to_output"
## Write output
write.table(resdata_KO1_padj, file.path(outDir, "J774_padj.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_KO1_padj_up_log2FC0, file.path(outDir, "J774_padj_up_log2FC0.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_KO1_padj_up_log2FC1, file.path(outDir, "J774_padj_up_log2FC1.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_KO1_padj_down_log2FC0, file.path(outDir, "J774_padj_down_log2FC0.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(resdata_KO1_padj_down_log2FC1, file.path(outDir, "J774_padj_down_log2FC1.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
```

## Filter out entries that have a padj < 0.5 (for creating list of constitutively expressed genes) and write
```{r}
resdata_KO1_insig <- subset(resdata_KO1, padj >= 0.75 & log2FoldChange > -0.5 & log2FoldChange < 0.5)
## Set output directory
outDir <- "/Users/%u/path_to_output"
## Write output
write.table(resdata_KO1_insig, file.path(outDir, "J774_insig_padj_0.75.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
```

# Make volcano plots
```{r}
library(ggplot2)
library(ggrepel)
library(dplyr)
padj <- res_KO1[,c(6)]
head(padj)
hist(padj)
hist(-log10(padj))
log2FoldChange <- res_KO1[, c(2)]
hist(log2FoldChange)
new <- res_KO1[, c(2,6)]
head(new)
new = as.data.frame(new)
value1 = subset(new, padj < 0.05 & abs(log2FoldChange) < 1)
value2 = subset(new, padj > 0.05 & abs(log2FoldChange) > 1)
value3 = subset(new, padj < 0.05 & abs(log2FoldChange) > 1)
dim(value3)
nrow(value3)
head(value3)
## these additional subsets are for the text labels 
value4 = subset(new, padj<0.0000000000000000000000000000005 & abs(log2FoldChange) > 4)
value4
dim(value4)
value5 = subset(new, padj<0.0000000000000000000000000000005 & log2FoldChange < -2.5)
value5
dim(value5)
## plot it
ggplot(new, aes(log2FoldChange, -log10(padj)), colour="grey") +
  scale_color_discrete(name = 'Labels') +
  theme_bw() + 
  ggtitle("Volcano plot: J774 mIFNG-treated WT (pooled) vs KO1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = "-log10 adjusted pvalue", x = "log2 fold change") +
  ## Set all dots color to grey
  geom_point(data=new, colour = "grey") + 
  ## If pvalue < 0.05, change dot color to green
  geom_point(data=new[which(new$padj<0.05),], colour = "springgreen2") + 
  ## If log2FC > 1, change dot color to orange
  geom_point(data=new[which(abs(new$log2FoldChange)>1),], colour = "darkgoldenrod1") +
  ## If both, change dot color to blue
  geom_point(data = new[which(abs(new$log2FoldChange)>1 & new$padj<0.05),], colour = "royalblue1") +
  # Add text label for interesting outliers
  geom_text_repel(data = value4, mapping = aes(log2FoldChange, -log10(padj), label = rownames(value4)), size = 3, force = 1, max.overlaps = Inf) +
  geom_text_repel(data = value5, mapping = aes(log2FoldChange, -log10(padj), label = rownames(value5)), size = 3, force = 1, max.overlaps = Inf)
dev.copy(pdf, height = 12, width = 12, pointsize=4, file.path(pngDir, "volcano.pdf"))
dev.off()
```
