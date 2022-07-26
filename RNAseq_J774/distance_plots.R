## Script for distance plots 

## Load libraries
library("dplyr")
library("ggplot2")
library("ggrepel")

## Set directories
## These should be the directories that host your count tables
inDir <- "/Users/%u/path_to_deseq2_output_table"
outDir <- "/Users/%u/path_to_deseq2_output_table/png"

## Read in the table
df <- read.table(file.path(inDir, "deseq2_table_including_geneCoords_chr12.txt"), header = TRUE, sep="\t")
head(df)

## Set first column (gene names) to be the rownames
rownames(df) <- df[,1]
head(df)

## Remove the first column (genename) and remote rows with "NA"
df <- df[ ,2:ncol(df)]
df <- na.omit(df)
head(df)

## Plot log2FC vs coordinates over entirety of Chr12
ggplot(df, aes(geneStart, log2FoldChange), color = "grey") +
  scale_color_discrete(name = "Labels") +
  theme_bw() +
  geom_hline(aes(yintercept = 0), size = 0.2) +
  ggtitle("Untreated WT vs mIFNG Treated WT J774s") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # change plot margins (order is top, right, bottom, left) and set font size
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(y = "Log2 Fold Change", x = "Genome Coordinates (bp)") +
  # If padj>0.05, color grey
  geom_point(data = df[which(df$padj >= 0.05), ], size = 1.5, colour = "grey", alpha = 0.7) +
  # add B2 annotation
  annotate("text", x = 115000000, y = 0.3, label = "B2 SINE", size = 4, fontface = 2, family="Helvetica") +
  annotate("rect", xmin=104742467-1500000, xmax = 104742467+1500000, ymin=-0.3, ymax=0.3, alpha=0.99, fill="black") +
  # If pvalue<0.05 & log2FC<0, increase dot size and make red
  geom_point(data = df[which(df$padj < 0.05 & df$log2FoldChange < 0), ], size = 3, colour = "red", alpha = 0.7) +
  # If pvalue<0.05 & log2FC>0, increase dot size and make blue 
  geom_point(data = df[which(df$padj < 0.05 & df$log2FoldChange > 0), ], size = 3, colour = "blue", alpha = 0.7) +
  # annotation for scale bar
  annotate("segment", x=10000, xend=10000+5000000, y=-11.3, yend=-11.3, fill="black") +
  annotate("text", x=10000+2500000, y=-10.6, label="5 Mb", size=3, family="Helvetica") +
  # Label genes within 1.5Mb of B2 SINE
  geom_label_repel(data = df[which(df$padj < 0.05 & df$geneStart > (104742467-1500000) & df$geneStart < (104742467+1500000)),],
                   mapping = aes(geneStart, 
                                 log2FoldChange,label = rownames(df[which(df$padj < 0.05 & df$geneStart > (104742467-1500000) & df$geneStart < (104742467+1500000)),])), 
                   size = 4, 
                   force = 10, 
                   max.overlaps = Inf, 
                   alpha = 0.8, ) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(expand = c(0,0), limits = c(-12,12)) +
  scale_size_continuous(range = c(0,7)) +
  theme(text=element_text(size=12, family="Helvetica"))

## Plot log2FC over 5Mb window centered on B2 SINE element on Chr12
ggplot(df, aes(geneStart, log2FoldChange), color = "grey") +
  scale_color_discrete(name = "Labels") +
  theme_bw() +
  geom_hline(aes(yintercept = 0), size = 0.2) +
  ggtitle("Untreated WT vs mIFNG Treated WT J774s") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # change plot margins (order is top, right, bottom, left) and set font size
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(y = "Log2 Fold Change", x = "Genome Coordinates (bp)") +
  # If padj>0.05, color grey
  geom_point(data = df[which(df$padj >= 0.05), ], size = 1.5, colour = "grey", alpha = 0.7) +
  # add B2 annotation
  annotate("text", x = 104742467+75000, y = -1, label = "B2 SINE", size = 4, fontface = 2, family="Helvetica") +
  annotate("rect", xmin=104742467-150000, xmax = 104742467+150000, ymin=-0.3, ymax=0.3, alpha=0.99, fill="black") +
  # If pvalue<0.05 & log2FC<0, increase dot size and make red 
  geom_point(data = df[which(df$padj < 0.05 & df$log2FoldChange < 0), ], size = 3, colour = "red", alpha = 0.7) +
  # If pvalue<0.05 & log2FC>0, increase dot size and make blue
  geom_point(data = df[which(df$padj < 0.05 & df$log2FoldChange > 0), ], size = 3, colour = "blue", alpha = 0.7) +
  # annotation for scale bar
  annotate("segment", x=104742467-2500000+100000, xend=104742467-2500000+100000+1500000, y=-11.5, yend=-11.5, colour="black") +
  annotate("text", x=104742467-2500000+100000+750000, y=-11, label="1.5 Mb", size=3, family="Helvetica") +
  # Label genes within 1.5Mb of B2 SINE
  geom_label_repel(data = df[which(df$padj < 0.05),],
                   mapping = aes(geneStart, 
                                 log2FoldChange,label = rownames(df[which(df$padj < 0.05),])), 
                   size = 4, 
                   force = 10, 
                   max.overlaps = Inf, 
                   alpha = 0.8, ) +
  scale_x_continuous(labels = scales::comma, limits = c(104742467-2500000,104742467+2500000)) +
  scale_y_continuous(expand = c(0,0), limits = c(-12,12)) +
  scale_size_continuous(range = c(0,7)) +
  theme(text=element_text(size=12, family="Helvetica"))
