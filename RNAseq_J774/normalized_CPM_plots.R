# R script for Serpina3g normalized counts from all samples DESeq Table
library(ggplot2)
library(ggpubr)

# Clear workspace
rm(list = ls())


# create an R table containing the conditions and normalized values for each condition
# can also import a table if you already have one made


condition <- factor(c(rep("UT WT", 3), rep("UT KO1", 3), rep("UT KO2", 3), rep("mIFNG WT", 3), rep("mIFNG KO1", 3), rep("mIFNG KO2", 3)))
condition <- data.frame(condition)
condition

# used normalized counts per million plus 1 and graphed via a log10 scale
countdata <- c(1+2.47724755765416, 1+1.22962072844653, 1+2.07693546916898, 1, 1, 1, 1, 1, 1, 1+3427.87638939083, 1+3693.05086795009, 1+3162.03857690326, 1+17.8891786639779, 1+16.5010575042858, 1+24.4792939539988, 1+10.0337412601264, 1+9.18390341168238, 1+5.34156106423287)

# bind the conditions to STS values
data <- cbind(countdata, condition)
data
colnames(data) <- c("value", "group")
data

data$group <- as.factor(data$group)

all_comparisons <- list( c("UT KO1","UT WT"), c("UT KO2","UT WT"), c("mIFNG KO1","mIFNG WT"), c("mIFNG KO2","mIFNG WT"), c("UT KO1","mIFNG KO1"), c("UT KO2","mIFNG KO2"), c("UT WT","mIFNG WT"))
all_comparisons

# make a boxplot with significance stars
ggplot(data, aes(x=group, y=value, fill=group)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0), size=3) +
  theme(plot.title = element_text(hjust = 0.1, size=12), axis.text=element_text(size=10)) +
  labs(y="Serpina3g Normalized Counts (CPM)", x = "group") +
  stat_compare_means(comparisons = all_comparisons, method = "t.test", label = "p.signif")

# calculate stats
library(Rmisc)
data_stats <- summarySE(data, measurevar="value", groupvars=c("group"))
data_stats

# reorder group categories so that dmso is first
data$group <- factor(data$group, levels = c('UT WT','mIFNG WT','UT KO1', 'mIFNG KO1', 'UT KO2', 'mIFNG KO2'))
data

# for both data and data_stats
data_stats$group <- factor(data_stats$group, levels = c('UT WT','mIFNG WT','UT KO1', 'mIFNG KO1', 'UT KO2', 'mIFNG KO2'))
data_stats

# Do a t-test to check it's significant
ggplot(data, aes(x=group, y=value)) +
  geom_bar(position=position_dodge(), stat="identity") +
  stat_compare_means(aes(x=group, y=value), comparisons = all_comparisons, method = "t.test", label = "p.signif") +
  theme_classic(base_size = 10) +
  geom_point(position="jitter", size=3) +
  scale_y_continuous(trans='log10', limits=c(1,1e7)) +
  labs(y="Serpina3g Normalized Counts (CPM)") +
  theme(axis.text.x = element_text(angle=-45, vjust=0.9, hjust=0))

# Add standard error bars to bar graph
ggplot(data_stats, aes(x=group, y=value)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2, position=position_dodge(.9)) +
  theme_classic(base_size = 10) +
  scale_y_continuous(trans='log10', limits=c(1,1e7)) +
  labs(y="Serpina3g Normalized Counts (CPM)") +
  theme(axis.text.x = element_text(angle=-45, vjust=0.9, hjust=0)) 

# export both plots (size 7 by 4)
# overlay in illustrator
