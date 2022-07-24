# Load required packages
library(ggplot2)
library(dplyr)

#### ISGs ONLY ####

# Read in tables
bound_Mm2_ISG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm2_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_nearest_2h_vs_UT_ISG_top750", header = FALSE, sep = "\t")
bound_Mm2_ISG <- bound_Mm2_ISG[,c(1:3, 13)]
colnames(bound_Mm2_ISG) <- c("chr", "start", "end", "distance")
head(bound_Mm2_ISG)

shuffle_Mm2_ISG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm2_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_shuffled_nearest_2h_vs_UT_ISG_top750", header = FALSE, sep = "\t")
shuffle_Mm2_ISG <- shuffle_Mm2_ISG[,c(1:3, 13)]
colnames(shuffle_Mm2_ISG) <- c("chr", "start", "end", "distance")
head(shuffle_Mm2_ISG)

unbound_Mm2_ISG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm2_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random2122_nearest_2h_vs_UT_ISG_top750", header = FALSE, sep = "\t")
unbound_Mm2_ISG <- unbound_Mm2_ISG[,c(1:3, 13)]
colnames(unbound_Mm2_ISG) <- c("chr", "start", "end", "distance")
head(unbound_Mm2_ISG)

unbound_Mm1a_ISG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm1a_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random2122_nearest_2h_vs_UT_ISG_top750", header = FALSE, sep = "\t")
unbound_Mm1a_ISG <- unbound_Mm1a_ISG[,c(1:3, 13)]
colnames(unbound_Mm1a_ISG) <- c("chr", "start", "end", "distance")
head(unbound_Mm1a_ISG)

unbound_Mm1t_ISG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm1t_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random2122_nearest_2h_vs_UT_ISG_top750", header = FALSE, sep = "\t")
unbound_Mm1t_ISG <- unbound_Mm1t_ISG[,c(1:3, 13)]
colnames(unbound_Mm1t_ISG) <- c("chr", "start", "end", "distance")
head(unbound_Mm1t_ISG)

# Subset
bound_Mm2_ISG_subset <- subset(bound_Mm2_ISG, distance < 500000)
head(bound_Mm2_ISG_subset)

shuffle_Mm2_ISG_subset <- subset(shuffle_Mm2_ISG, distance < 500000)
head(shuffle_Mm2_ISG_subset)

unbound_Mm2_ISG_subset <- subset(unbound_Mm2_ISG, distance < 500000)
head(unbound_Mm2_ISG_subset)

unbound_Mm1a_ISG_subset <- subset(unbound_Mm1a_ISG, distance < 500000)
head(unbound_Mm1a_ISG_subset)

unbound_Mm1t_ISG_subset <- subset(unbound_Mm1t_ISG, distance < 500000)
head(unbound_Mm1t_ISG_subset)

# Convert from bp to kb & add group information
bound_Mm2_ISG_distanceOnly <- bound_Mm2_ISG_subset$distance / 1000
bound_Mm2_ISG_distanceOnly <- data.frame(bound_Mm2_ISG_distanceOnly, rep("bound_B2_Mm2_ISG", nrow(bound_Mm2_ISG_subset)))
colnames(bound_Mm2_ISG_distanceOnly) <- c("distance", "group")

shuffle_Mm2_ISG_distanceOnly <- shuffle_Mm2_ISG_subset$distance / 1000
shuffle_Mm2_ISG_distanceOnly <- data.frame(shuffle_Mm2_ISG_distanceOnly, rep("shuffled_bound_B2_Mm2_ISG", nrow(shuffle_Mm2_ISG_subset)))
colnames(shuffle_Mm2_ISG_distanceOnly) <- c("distance", "group")

unbound_Mm2_ISG_distanceOnly <- unbound_Mm2_ISG_subset$distance / 1000
unbound_Mm2_ISG_distanceOnly <- data.frame(unbound_Mm2_ISG_distanceOnly, rep("unbound_B2_Mm2_ISG", nrow(unbound_Mm2_ISG_subset)))
colnames(unbound_Mm2_ISG_distanceOnly) <- c("distance", "group")

unbound_Mm1a_ISG_distanceOnly <- unbound_Mm1a_ISG_subset$distance / 1000
unbound_Mm1a_ISG_distanceOnly <- data.frame(unbound_Mm1a_ISG_distanceOnly, rep("unbound_B2_Mm1a_ISG", nrow(unbound_Mm1a_ISG_subset)))
colnames(unbound_Mm1a_ISG_distanceOnly) <- c("distance", "group")

unbound_Mm1t_ISG_distanceOnly <- unbound_Mm1t_ISG_subset$distance / 1000
unbound_Mm1t_ISG_distanceOnly <- data.frame(unbound_Mm1t_ISG_distanceOnly, rep("unbound_B2_Mm1t_ISG", nrow(unbound_Mm1t_ISG_subset)))
colnames(unbound_Mm1t_ISG_distanceOnly) <- c("distance", "group")

# Combine
combined_ISG <- rbind(bound_Mm2_ISG_distanceOnly, shuffle_Mm2_ISG_distanceOnly, unbound_Mm2_ISG_distanceOnly, unbound_Mm1a_ISG_distanceOnly, unbound_Mm1t_ISG_distanceOnly)

# Plot number of peaks (y) vs. distance from nearest TSS (x) (legend)
pdf(file = "nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_STAT1_nearest_2h_vs_UT_ISG_top750_freqHist.pdf", width = 7, height = 7)

ggplot(data = combined_IRG) +
  geom_freqpoly(binwidth = 20, 
                size = 1,
                aes(x = distance,
                    color = group,
                    y = ..count..),
                na.rm = TRUE) +
  labs(x = "Distance (kb) from nearest mouse platanitis BMDM ISG", y = "Number of B2-derived STAT1 platanitis BMDM peaks") +
  scale_color_manual(values = c("royalblue3", "black", "#A9A9A9", "#E5E4E2","#7393B3")) +
  theme_bw() +
  xlim(0, 500) +
  ylim(0, 100) +
  theme(axis.title.x = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        legend.title = element_blank()) 

dev.off()

#### IRGs ONLY ####

# Read in tables
bound_Mm2_IRG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm2_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_nearest_2h_vs_UT_IRG_top750", header = FALSE, sep = "\t")
bound_Mm2_IRG <- bound_Mm2_IRG[,c(1:3, 13)]
colnames(bound_Mm2_IRG) <- c("chr", "start", "end", "distance")
head(bound_Mm2_IRG)

shuffle_Mm2_IRG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm2_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_shuffled_nearest_2h_vs_UT_IRG_top750", header = FALSE, sep = "\t")
shuffle_Mm2_IRG <- shuffle_Mm2_IRG[,c(1:3, 13)]
colnames(shuffle_Mm2_IRG) <- c("chr", "start", "end", "distance")
head(shuffle_Mm2_IRG)

unbound_Mm2_IRG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm2_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random2122_nearest_2h_vs_UT_IRG_top750", header = FALSE, sep = "\t")
unbound_Mm2_IRG <- unbound_Mm2_IRG[,c(1:3, 13)]
colnames(unbound_Mm2_IRG) <- c("chr", "start", "end", "distance")
head(unbound_Mm2_IRG)

unbound_Mm1a_IRG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm1a_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random2122_nearest_2h_vs_UT_IRG_top750", header = FALSE, sep = "\t")
unbound_Mm1a_IRG <- unbound_Mm1a_IRG[,c(1:3, 13)]
colnames(unbound_Mm1a_IRG) <- c("chr", "start", "end", "distance")
head(unbound_Mm1a_IRG)

unbound_Mm1t_IRG <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm1t_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random2122_nearest_2h_vs_UT_IRG_top750", header = FALSE, sep = "\t")
unbound_Mm1t_IRG <- unbound_Mm1t_IRG[,c(1:3, 13)]
colnames(unbound_Mm1t_IRG) <- c("chr", "start", "end", "distance")
head(unbound_Mm1t_IRG)

# Subset
bound_Mm2_IRG_subset <- subset(bound_Mm2_IRG, distance < 500000)
head(bound_Mm2_IRG_subset)

shuffle_Mm2_IRG_subset <- subset(shuffle_Mm2_IRG, distance < 500000)
head(shuffle_Mm2_IRG_subset)

unbound_Mm2_IRG_subset <- subset(unbound_Mm2_IRG, distance < 500000)
head(unbound_Mm2_IRG_subset)

unbound_Mm1a_IRG_subset <- subset(unbound_Mm1a_IRG, distance < 500000)
head(unbound_Mm1a_IRG_subset)

unbound_Mm1t_IRG_subset <- subset(unbound_Mm1t_IRG, distance < 500000)
head(unbound_Mm1t_IRG_subset)

# Convert from bp to kb & add group information
bound_Mm2_IRG_distanceOnly <- bound_Mm2_IRG_subset$distance / 1000
bound_Mm2_IRG_distanceOnly <- data.frame(bound_Mm2_IRG_distanceOnly, rep("bound_B2_Mm2_IRG", nrow(bound_Mm2_IRG_subset)))
colnames(bound_Mm2_IRG_distanceOnly) <- c("distance", "group")

shuffle_Mm2_IRG_distanceOnly <- shuffle_Mm2_IRG_subset$distance / 1000
shuffle_Mm2_IRG_distanceOnly <- data.frame(shuffle_Mm2_IRG_distanceOnly, rep("shuffled_bound_B2_Mm2_IRG", nrow(shuffle_Mm2_IRG_subset)))
colnames(shuffle_Mm2_IRG_distanceOnly) <- c("distance", "group")

unbound_Mm2_IRG_distanceOnly <- unbound_Mm2_IRG_subset$distance / 1000
unbound_Mm2_IRG_distanceOnly <- data.frame(unbound_Mm2_IRG_distanceOnly, rep("unbound_B2_Mm2_IRG", nrow(unbound_Mm2_IRG_subset)))
colnames(unbound_Mm2_IRG_distanceOnly) <- c("distance", "group")

unbound_Mm1a_IRG_distanceOnly <- unbound_Mm1a_IRG_subset$distance / 1000
unbound_Mm1a_IRG_distanceOnly <- data.frame(unbound_Mm1a_IRG_distanceOnly, rep("unbound_B2_Mm1a_IRG", nrow(unbound_Mm1a_IRG_subset)))
colnames(unbound_Mm1a_IRG_distanceOnly) <- c("distance", "group")

unbound_Mm1t_IRG_distanceOnly <- unbound_Mm1t_IRG_subset$distance / 1000
unbound_Mm1t_IRG_distanceOnly <- data.frame(unbound_Mm1t_IRG_distanceOnly, rep("unbound_B2_Mm1t_IRG", nrow(unbound_Mm1t_IRG_subset)))
colnames(unbound_Mm1t_IRG_distanceOnly) <- c("distance", "group")

# Combine
combined_IRG <- rbind(bound_Mm2_IRG_distanceOnly, shuffle_Mm2_IRG_distanceOnly, unbound_Mm2_IRG_distanceOnly, unbound_Mm1a_IRG_distanceOnly, unbound_Mm1t_IRG_distanceOnly)

# Plot number of peaks (y) vs. distance from nearest TSS (x) (legend)
pdf(file = "nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_STAT1_nearest_2h_vs_UT_IRG_top750_freqHist.pdf", width = 7, height = 7)

ggplot(data = combined_IRG) +
  geom_freqpoly(binwidth = 20, 
                size = 1,
                aes(x = distance,
                    color = group,
                    y = ..count..),
                na.rm = TRUE) +
  labs(x = "Distance (kb) from nearest mouse platanitis BMDM IRG", y = "Number of B2-derived STAT1 platanitis BMDM peaks") +
  scale_color_manual(values = c("royalblue3", "black", "#A9A9A9", "#E5E4E2","#7393B3")) +
  theme_bw() +
  xlim(0, 500) +
  ylim(0, 100) +
  theme(axis.title.x = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        legend.title = element_blank()) 

dev.off()

#### NONRESPONSIVE GENES ONLY ####

# Read in tables
bound_Mm2_nonresponsive <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm2_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_nearest_2h_vs_UT_nonresponsive_random750", header = FALSE, sep = "\t")
bound_Mm2_nonresponsive <- bound_Mm2_nonresponsive[,c(1:3, 13)]
colnames(bound_Mm2_nonresponsive) <- c("chr", "start", "end", "distance")
head(bound_Mm2_nonresponsive)

shuffle_Mm2_nonresponsive <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm2_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_shuffled_nearest_2h_vs_UT_nonresponsive_random750", header = FALSE, sep = "\t")
shuffle_Mm2_nonresponsive <- shuffle_Mm2_nonresponsive[,c(1:3, 13)]
colnames(shuffle_Mm2_nonresponsive) <- c("chr", "start", "end", "distance")
head(shuffle_Mm2_nonresponsive)

unbound_Mm2_nonresponsive <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm2_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random2122_nearest_2h_vs_UT_nonresponsive_random750", header = FALSE, sep = "\t")
unbound_Mm2_nonresponsive <- unbound_Mm2_nonresponsive[,c(1:3, 13)]
colnames(unbound_Mm2_nonresponsive) <- c("chr", "start", "end", "distance")
head(unbound_Mm2_nonresponsive)

unbound_Mm1a_nonresponsive <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm1a_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random2122_nearest_2h_vs_UT_nonresponsive_random750", header = FALSE, sep = "\t")
unbound_Mm1a_nonresponsive <- unbound_Mm1a_nonresponsive[,c(1:3, 13)]
colnames(unbound_Mm1a_nonresponsive) <- c("chr", "start", "end", "distance")
head(unbound_Mm1a_nonresponsive)

unbound_Mm1t_nonresponsive <- read.table("nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_Mm1t_null_plat_BMDM_WT_IFNG_1.5h_STAT1_intersected_random2122_nearest_2h_vs_UT_nonresponsive_random750", header = FALSE, sep = "\t")
unbound_Mm1t_nonresponsive <- unbound_Mm1t_nonresponsive[,c(1:3, 13)]
colnames(unbound_Mm1t_nonresponsive) <- c("chr", "start", "end", "distance")
head(unbound_Mm1t_nonresponsive)

# Subset
bound_Mm2_nonresponsive_subset <- subset(bound_Mm2_nonresponsive, distance < 500000)
head(bound_Mm2_nonresponsive_subset)

shuffle_Mm2_nonresponsive_subset <- subset(shuffle_Mm2_nonresponsive, distance < 500000)
head(shuffle_Mm2_nonresponsive_subset)

unbound_Mm2_nonresponsive_subset <- subset(unbound_Mm2_nonresponsive, distance < 500000)
head(unbound_Mm2_nonresponsive_subset)

unbound_Mm1a_nonresponsive_subset <- subset(unbound_Mm1a_nonresponsive, distance < 500000)
head(unbound_Mm1a_nonresponsive_subset)

unbound_Mm1t_nonresponsive_subset <- subset(unbound_Mm1t_nonresponsive, distance < 500000)
head(unbound_Mm1t_nonresponsive_subset)

# Convert from bp to kb & add group information
bound_Mm2_nonresponsive_distanceOnly <- bound_Mm2_nonresponsive_subset$distance / 1000
bound_Mm2_nonresponsive_distanceOnly <- data.frame(bound_Mm2_nonresponsive_distanceOnly, rep("bound_B2_Mm2_nonresponsive", nrow(bound_Mm2_nonresponsive_subset)))
colnames(bound_Mm2_nonresponsive_distanceOnly) <- c("distance", "group")

shuffle_Mm2_nonresponsive_distanceOnly <- shuffle_Mm2_nonresponsive_subset$distance / 1000
shuffle_Mm2_nonresponsive_distanceOnly <- data.frame(shuffle_Mm2_nonresponsive_distanceOnly, rep("shuffled_bound_B2_Mm2_nonresponsive", nrow(shuffle_Mm2_nonresponsive_subset)))
colnames(shuffle_Mm2_nonresponsive_distanceOnly) <- c("distance", "group")

unbound_Mm2_nonresponsive_distanceOnly <- unbound_Mm2_nonresponsive_subset$distance / 1000
unbound_Mm2_nonresponsive_distanceOnly <- data.frame(unbound_Mm2_nonresponsive_distanceOnly, rep("unbound_B2_Mm2_nonresponsive", nrow(unbound_Mm2_nonresponsive_subset)))
colnames(unbound_Mm2_nonresponsive_distanceOnly) <- c("distance", "group")

unbound_Mm1a_nonresponsive_distanceOnly <- unbound_Mm1a_nonresponsive_subset$distance / 1000
unbound_Mm1a_nonresponsive_distanceOnly <- data.frame(unbound_Mm1a_nonresponsive_distanceOnly, rep("unbound_B2_Mm1a_nonresponsive", nrow(unbound_Mm1a_nonresponsive_subset)))
colnames(unbound_Mm1a_nonresponsive_distanceOnly) <- c("distance", "group")

unbound_Mm1t_nonresponsive_distanceOnly <- unbound_Mm1t_nonresponsive_subset$distance / 1000
unbound_Mm1t_nonresponsive_distanceOnly <- data.frame(unbound_Mm1t_nonresponsive_distanceOnly, rep("unbound_B2_Mm1t_nonresponsive", nrow(unbound_Mm1t_nonresponsive_subset)))
colnames(unbound_Mm1t_nonresponsive_distanceOnly) <- c("distance", "group")

# Combine
combined_nonresponsive <- rbind(bound_Mm2_nonresponsive_distanceOnly, shuffle_Mm2_nonresponsive_distanceOnly, unbound_Mm2_nonresponsive_distanceOnly, unbound_Mm1a_nonresponsive_distanceOnly, unbound_Mm1t_nonresponsive_distanceOnly)

# Plot number of peaks (y) vs. distance from nearest TSS (x) (legend)
pdf(file = "nearest_neighbor/platanitis_RNAseq_IFNG_2h_vs_UT/B2_STAT1_nearest_2h_vs_UT_nonresponsive_random750_freqHist.pdf", width = 7, height = 7)

ggplot(data = combined_nonresponsive) +
  geom_freqpoly(binwidth = 20, 
                size = 1,
                aes(x = distance,
                    color = group,
                    y = ..count..),
                na.rm = TRUE) +
  labs(x = "Distance (kb) from nearest mouse platanitis BMDM nonresponsive", y = "Number of B2-derived STAT1 platanitis BMDM peaks") +
  scale_color_manual(values = c("royalblue3", "black", "#A9A9A9", "#E5E4E2","#7393B3")) +
  theme_bw() +
  xlim(0, 500) +
  ylim(0, 100) +
  theme(axis.title.x = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, hjust = 0.5, vjust = 0.5, color = "black"),
        legend.title = element_blank()) 

dev.off()