# Load packages
library("dplyr", "tidyr", "plotrix", "ggplot2")

# Set input directory
# This should be the directory that hosts your TEtranscripts DESeq2 results tables
workingDir <- "/Users/path/to/working/directory"

# Read in TEtranscripts DESeq2 results
data_picc_2h <- read.csv(file.path(workingDir, "picc_BMDM_IFNG_2h_vs_UT_TEtranscripts.txt"), sep = "", header = FALSE)
data_plat_2h <- read.csv(file.path(workingDir, "picc_BMDM_IFNG_4h_vs_UT_TEtranscripts.txt"), sep = "", header = FALSE)
data_picc_4h <- read.csv(file.path(workingDir, "plat_BMDM_IFNG_2h_vs_UT_TEtranscripts.txt"), sep = "", header = FALSE)

# Filter all columns except those containing family names or counts
data_picc_2h <- select(data_picc_2h, !contains(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")))
data_picc_4h <- select(data_picc_4h, !contains(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")))
data_plat_2h <- select(data_plat_2h, !contains(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")))

# Rename sample columns to add dataset information
colnames(data_picc_2h) <- c("family", "picc_UT_R1", "picc_UT_R3", "picc_IFNG_2h_R1", "picc_IFNG_2h_R2", "picc_IFNG_2h_R3")
colnames(data_picc_4h) <- c("family", "picc_UT_R1", "picc_UT_R3", "picc_IFNG_4h_R1", "picc_IFNG_4h_R2", "picc_IFNG_4h_R3")
colnames(data_plat_2h) <- c("family", "plat_UT_R1", "plat_UT_R2", "plat_IFNG_2h_R1", "plat_IFNG_2h_R2", "plat_IFNG_2h_R3")

# Combine Piccolo et al. 2017 tables & remove redundant columns
data_picc_combined <- left_join(data_picc_2h, data_picc_4h, by = "family") %>% select(!contains(c(".y")))
colnames(data_picc_combined) <-  c("family", "picc_UT_R1", "picc_UT_R3", 
                                   "picc_IFNG_2h_R1", "picc_IFNG_2h_R2", "picc_IFNG_2h_R3", 
                                   "picc_IFNG_4h_R1", "picc_IFNG_4h_R2", "picc_IFNG_4h_R3")

# Convert from wide to long
data_picc_combined <- gather(data_picc_combined, sample, counts, picc_UT_R1:picc_IFNG_4h_R3, factor_key=TRUE)
data_plat_2h <- gather(data_plat_2h, sample, counts, plat_UT_R1:plat_IFNG_2h_R3, factor_key=TRUE)

# Add group column
data_picc_combined <- mutate(data_picc_combined, group = case_when(grepl("picc_UT", sample) ~ "picc_UT",
                                                                   grepl("picc_IFNG_2h", sample) ~ "picc_IFNG_2h",
                                                                   grepl("picc_IFNG_4h", sample) ~ "picc_IFNG_4h"))
data_plat_2h <- mutate(data_plat_2h, group = case_when(grepl("plat_UT", sample) ~ "plat_UT",
                                                       grepl("plat_IFNG_2h", sample) ~ "plat_IFNG_2h"))

# Define factor levels
data_picc_combined$group <- factor(data_picc_combined$group, levels = c("picc_UT", "picc_IFNG_2h", "picc_IFNG_4h"))
data_plat_2h$group <- factor(data_plat_2h$group, levels = c("plat_UT", "plat_IFNG_2h"))

# Subset for unbound B2_Mm2
data_picc_combined_B2_Mm2 <- filter(data_picc_combined, family == "Unbound_B2_Mm2")
data_plat_2h_B2_Mm2 <- filter(data_plat_2h, family == "Unbound_B2_Mm2")

# Subset for STAT1-bound B2_Mm2
data_picc_combined_STAT1_B2_Mm2 <- filter(data_picc_combined, family == "STAT1_B2_Mm2")
data_plat_2h_STAT1_B2_Mm2 <- filter(data_plat_2h, family == "STAT1_B2_Mm2")

# Get average for each treatment condition per group for unbound B2_Mm2
data_picc_combined_B2_Mm2 <- data_picc_combined_B2_Mm2 %>% group_by(group) %>% mutate(mean = mean(counts), se = std.error(counts)) %>% as.data.frame()
data_plat_2h_B2_Mm2 <- data_plat_2h_B2_Mm2 %>% group_by(group) %>% mutate(mean = mean(counts), se = std.error(counts)) %>% as.data.frame()

# Get average for each treatment condition per group for STAT1-bound B2_Mm2
data_picc_combined_STAT1_B2_Mm2 <- data_picc_combined_STAT1_B2_Mm2 %>% group_by(group) %>% mutate(mean = mean(counts), se = std.error(counts)) %>% as.data.frame()
data_plat_2h_STAT1_B2_Mm2 <- data_plat_2h_STAT1_B2_Mm2 %>% group_by(group) %>% mutate(mean = mean(counts), se = std.error(counts)) %>% as.data.frame()

# Import font
font_add("Arial", "/scratch/Users/coke6162/tmp/Fonts/Arial.ttf")
showtext_auto()

# Plot unbound B2_Mm2 (Piccolo et al. 2017)
picc_combined_B2_Mm2 <- ggplot(data = data_picc_combined_B2_Mm2, 
                               aes(x = group, y = mean, fill = group)) +
  geom_bar(position = position_dodge(0.1), stat = 'identity', width = 0.5) +
  geom_errorbar(mapping = aes(x = group, ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_point(aes(y = counts), position = position_dodge2(0.25, preserve = "total"), size = 1.5, shape = 16) +
  labs(x = NULL, y = "Normalized counts") +
  scale_fill_manual(values = c("#C1C0C0", "#9b9ecb", "#3C55A4")) +
  theme_classic(base_size = 10, base_family = "Arial") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 22500)) +
  theme(axis.title.y = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        legend.position = "none")
picc_combined_B2_Mm2
ggsave("bar_picc_2h_and_4h_vs_UT_Unbound_B2_Mm2.pdf", plot = picc_combined_B2_Mm2, height = 2, width = 2.75, units = "in", useDingbats = FALSE)

# Plot unbound B2_Mm2 (Platanitis et al. 2019)
plat_2h_B2_Mm2 <- ggplot(data = data_plat_2h_B2_Mm2, 
               aes(x = group, y = mean, fill = group)) +
  geom_bar(position = position_dodge(0.1), stat = 'identity', width = 0.5) +
  geom_errorbar(mapping = aes(x = group, ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_point(aes(y = counts), position = position_dodge2(0.25, preserve = "total"), size = 1.5, shape = 16) +
  labs(x = NULL, y = "Normalized counts") +
  scale_fill_manual(values = c("#C1C0C0", "#9b9ecb")) +
  theme_classic(base_size = 10, base_family = "Arial") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 8000)) +
  theme(axis.title.y = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        legend.position = "none")
plat_2h_B2_Mm2
ggsave("bar_plat_2h_vs_UT_Unbound_B2_Mm2.pdf", plot = plat_2h_B2_Mm2, height = 2, width = 2.75, units = "in", useDingbats = FALSE)

# Plot STAT1-bound B2_Mm2 (Piccolo et al. 2017)
picc_combined_STAT1_B2_Mm2 <- ggplot(data = data_picc_combined_STAT1_B2_Mm2, 
                                     aes(x = group, y = mean, fill = group)) +
  geom_bar(position = position_dodge(0.1), stat = 'identity', width = 0.5) +
  geom_errorbar(mapping = aes(x = group, ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_point(aes(y = counts), position = position_dodge2(0.25, preserve = "total"), size = 1.5, shape = 16) +
  labs(x = NULL, y = "Normalized counts") +
  scale_fill_manual(values = c("#C1C0C0", "#9b9ecb", "#3C55A4")) +
  theme_classic(base_size = 10, base_family = "Arial") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 4250)) +
  theme(axis.title.y = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        legend.position = "none")
picc_combined_STAT1_B2_Mm2
ggsave("bar_picc_2h_and_4h_vs_UT_STAT1_B2_Mm2.pdf", plot = picc_combined_STAT1_B2_Mm2, height = 2, width = 2.75, units = "in", useDingbats = FALSE)

# Plot STAT1-bound B2_Mm2 (Platanitis et al. 2019)
plat_2h_STAT1_B2_Mm2 <- ggplot(data = data_plat_2h_STAT1_B2_Mm2, 
                                     aes(x = group, y = mean, fill = group)) +
  geom_bar(position = position_dodge(0.1), stat = 'identity', width = 0.5) +
  geom_errorbar(mapping = aes(x = group, ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_point(aes(y = counts), position = position_dodge2(0.25, preserve = "total"), size = 1.5, shape = 16) +
  labs(x = NULL, y = "Normalized counts") +
  scale_fill_manual(values = c("#C1C0C0", "#9b9ecb")) +
  theme_classic(base_size = 10, base_family = "Arial") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 250)) +
  theme(axis.title.y = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        legend.position = "none")
plat_2h_STAT1_B2_Mm2
ggsave("bar_plat_2h_vs_UT_STAT1_B2_Mm2.pdf", plot = plat_2h_STAT1_B2_Mm2, height = 2, width = 2.75, units = "in", useDingbats = FALSE)