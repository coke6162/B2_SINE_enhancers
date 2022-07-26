# This is an example script of what was used in R Studio to create the differential expression qPCR plots

library(ggplot2)
library(ggpubr)

# Clear workspace
rm(list = ls())


# create an R table containing the conditions and normalized values for each condition
# can also import a table if you already have one made


condition <- factor(c(rep("UT WT1",3), rep("mIFNG WT1",3), rep("UT WT2",3), rep("mIFNG WT2",3), rep("UT WT3",3), rep("mIFNG WT3",3), rep("UT KO1",3), rep("mIFNG KO1",3), rep("UT KO2",3), rep("mIFNG KO2",3)))
condition <- data.frame(condition)
condition

countdata <- c(0.945204754, 1.028303107, 0.899408047, 1.403626906, 1.204685099, 1.416580769, 0.852727738, 0.69409343, 0.695597853, 1.410968086, 1.672755674, 1.505040265, 0.850973528, 1.167002564, 1.290197241,
               1.692917679, 1.506396851, 1.545405804, 0.678882001, 0.717726928, 0.630145209, 0.776477305, 0.724617612, 0.531877304, 0.812664702, 0.820105545, 0.868454252, 0.814679607, 0.947914034, 1.020852177)

# bind the conditions to STS values
data <- cbind(countdata, condition)
data
colnames(data) <- c("value", "group")
data

data$group <- as.factor(data$group)

all_comparisons <- list( c("UT WT1","mIFNG WT1"), c("UT WT2","mIFNG WT2"), c("UT WT3","mIFNG WT3"), c("UT KO1","mIFNG KO1"), c("UT KO2","mIFNG KO2"), c("mIFNG KO1","mIFNG WT1"), c("mIFNG KO1","mIFNG WT2"), c("mIFNG KO1","mIFNG WT3"), c("mIFNG KO2", "mIFNG WT1"), c("mIFNG KO2", "mIFNG WT2"), c("mIFNG KO2", "mIFNG WT3"))
all_comparisons

# make a boxplot with significance stars
ggplot(data, aes(x=group, y=value, fill=group)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0), size=3) +
  theme(plot.title = element_text(hjust = 0.1, size=12), axis.text=element_text(size=10)) +
  theme(axis.text.x = element_text(angle=-45, vjust=0.9, hjust=0.2)) +
  labs(y="Dicer1 DeltaCq Expression (Normalized to CTCF)", x = "group") +
  stat_compare_means(comparisons = all_comparisons, method = "t.test", label = "p.signif")

# calculate stats
library(Rmisc)
data_stats <- summarySE(data, measurevar="value", groupvars=c("group"))
data_stats

# reorder group categories so that dmso is first
data$group <- factor(data$group, levels = c('UT WT1','mIFNG WT1', 'UT WT2','mIFNG WT2', 'UT WT3','mIFNG WT3', 'UT KO1', 'mIFNG KO1', 'UT KO2', 'mIFNG KO2'))
data

# for both data and data_stats
data_stats$group <- factor(data_stats$group, levels = c('UT WT1','mIFNG WT1', 'UT WT2','mIFNG WT2', 'UT WT3','mIFNG WT3', 'UT KO1', 'mIFNG KO1', 'UT KO2', 'mIFNG KO2'))
data_stats

# Do a t-test to check it's significant
ggplot(data, aes(x=group, y=value)) +
  geom_bar(position=position_dodge(), stat="identity") +
  stat_compare_means(aes(x=group, y=value), comparisons = all_comparisons, method = "t.test", var.equal = "TRUE", label = "p.signif") +
  theme_classic(base_size = 10) +
  geom_point(position="jitter", size=3) +
  labs(y="Dicer1 DeltaCq Expression (Normalized to CTCF)") +
  theme(axis.text.x = element_text(angle=-45, vjust=0.9, hjust=0.2)) +
  ylim(0,3.5)

# Add standard error bars to bar graph
ggplot(data_stats, aes(x=group, y=value)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2, position=position_dodge(.9)) +
  theme_classic(base_size = 10) +
  labs(y="Dicer1 DeltaCq Expression (Normalized to CTCF)") +
  theme(axis.text.x = element_text(angle=-45, vjust=0.9, hjust=0.2)) +
  ylim(0,3.5)
