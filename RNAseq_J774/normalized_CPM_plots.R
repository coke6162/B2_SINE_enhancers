library(ggplot2)
library(ggpubr)

# Clear workspace
rm(list = ls())


# create an R table containing the conditions and normalized values for each condition
# can also import a table if you already have one made


condition <- factor(c(rep("UT WT1", 3), rep("UT KO1", 3), rep("UT KO2", 3), rep("mIFNG WT1", 3), rep("mIFNG KO1", 3), rep("mIFNG KO2", 3), rep("UT WT2", 3), rep("UT WT3", 3), rep("mIFNG WT2", 3), rep("mIFNG WT3", 3)))
condition <- data.frame(condition)
condition

countdata <- c(2205.73759948494, 2515.51929803945, 2484.40770629049, 1609.96289057184, 1658.33046789426, 1600.15865496344, 2485.23061938634, 2368.40056105501, 2306.69307171249, 3169.81480306144, 2988.27161873072, 3021.7688020462, 1842.80938723078, 1628.10667821876, 1659.36883935643, 2205.88788755095, 2377.20656872349, 2211.44821893401, 2853.29, 2962.83, 2947.31, 2855.64, 2861.43, 2755.74, 3634.06, 3591.447299, 3659.34, 3472.57,3930.15,3700.77)

# bind the conditions to STS values
data <- cbind(countdata, condition)
data
colnames(data) <- c("value", "group")
data

data$group <- as.factor(data$group)

all_comparisons <- list( c("UT KO1","UT WT1"), c("UT KO2","UT WT1"), c("mIFNG KO1","mIFNG WT1"), c("mIFNG KO2","mIFNG WT1"), c("UT KO1","mIFNG KO1"), c("UT KO2","mIFNG KO2"), c("UT WT1","mIFNG WT1"), c("UT WT2","mIFNG WT2"), c("UT WT3","mIFNG WT3"), c("UT KO1","UT WT2"), c("UT KO2","UT WT2"), c("UT KO1","UT WT3"), c("UT KO2","UT WT3"), c("mIFNG KO1","mIFNG WT2"), c("mIFNG KO2","mIFNG WT2"), c("mIFNG KO1","mIFNG WT3"), c("mIFNG KO2","mIFNG WT3"))
all_comparisons

# make a boxplot with significance stars
ggplot(data, aes(x=group, y=value, fill=group)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0), size=3) +
  theme(plot.title = element_text(hjust = 0.1, size=12), axis.text=element_text(size=10)) +
  labs(y="Dicer1 Normalized Counts", x = "group") +
  stat_compare_means(comparisons = all_comparisons, method = "t.test", label = "p.signif")

# calculate stats
library(Rmisc)
data_stats <- summarySE(data, measurevar="value", groupvars=c("group"))
data_stats

# reorder group categories so that dmso is first
data$group <- factor(data$group, levels = c('UT WT1','mIFNG WT1','UT WT2','mIFNG WT2','UT WT3','mIFNG WT3','UT KO1', 'mIFNG KO1', 'UT KO2', 'mIFNG KO2'))
data

# for both data and data_stats
data_stats$group <- factor(data_stats$group, levels = c('UT WT1','mIFNG WT1','UT WT2','mIFNG WT2','UT WT3','mIFNG WT3','UT KO1', 'mIFNG KO1', 'UT KO2', 'mIFNG KO2'))
data_stats

# Do a t-test to check it's significant
ggplot(data, aes(x=group, y=value)) +
  geom_bar(position=position_dodge(), stat="identity") +
  stat_compare_means(aes(x=group, y=value), comparisons = all_comparisons, method = "t.test", label = "p.signif") +
  theme_classic(base_size = 10) +
  geom_point(position="jitter", size=3) +
  labs(y="Dicer1 Normalized Counts") +
  theme(axis.text.x = element_text(angle=-45, vjust=0.9, hjust=0)) +
  ylim(0,4500)

# Add standard error bars to bar graph
ggplot(data_stats, aes(x=group, y=value)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2, position=position_dodge(.9)) +
  theme_classic(base_size = 10) +
  labs(y="Dicer1 Normalized Counts") +
  theme(axis.text.x = element_text(angle=-45, vjust=0.9, hjust=0)) +
  ylim(0,4500)

# export both plots (size 7 by 4)
# overlay in illustrator