# Load required packages
# Credit for function: https://gist.github.com/smithdanielle/9913897
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("dplyr", "ggplot2", "reshape2", "cowplot", "scales")
check.packages(packages)

# Read in input file & rename columns
data <- read.table("box_and_whisker/B2_GAS_pval_all_summarized.txt", sep = "\t", header = FALSE)
colnames(data) <- c("group", "pval")

# Calculate sample size for each group
data_n <- data %>%
  group_by(group) %>%
  mutate(n = n()) %>%
  mutate(label = paste0(group,'\nn = ',n))

# Box and whisker
whisker <- ggplot(data_n, aes(x = reorder(factor(label), pval, FUN = median),
                   y = -log10(pval),
                   color = "black",
                   fill = label)) +
  geom_boxplot(width = 5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.44), size = 0.05, alpha = 0.1) +
  theme_bw() +
  scale_color_manual(values=c("black"), guide = "none") +
  scale_fill_manual(values=c("#d7191c", "#fdae61", "#e6e6ac", "#abdda4", "#2b83ba"), guide = "none") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 6.5)) +
  ylab("-Log10 FIMO p-value") +
  theme(axis.title.y = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none")

# Save
options(bitmapType = 'cairo')
ggsave("box_and_whisker/B2_GAS_whisker.pdf",
       plot = whisker,
       device = png(),
       width = unit(4, "in"),
       height = unit(6, "in"),
       limitsize = FALSE)