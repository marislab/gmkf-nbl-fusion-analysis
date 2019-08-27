# Author: Komal S. Rathi
# Date: 08/19/2019
# Function: Create Summary Plots
# Step 4
source('~/Projects/DGD_Mendelian_RNASeq/R/pubTheme.R')

# Line plot
df <- data.frame(Filter = c("F0","F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8"), 
                 Label = c("Combined Calls", "Filter Circular RNAs", 
                           "FusionAnnotator", "Filter Read-throughs", 
                           "Filter Duplicates", "Filter by Gene Biotype", 
                           "Annotation filters", "Expression Filter", "Outlier Analysis"),
                 Val = c(5062, 4846, 4716, 3609, 3596, 3386, 2741, 2739, 364))
df$Sample <- "GMKF_NBL"


p <- ggplot(df, aes(x = Filter, y = Val, group = Sample)) + 
  geom_point(aes(color = Filter), pch = 18, size = 6) + 
  geom_line(linetype = "dotted") + 
  geom_text(aes(label = Val, hjust = 0.5, vjust = -1)) + theme_Publication3() + 
  theme(legend.position = "right") +
  scale_color_discrete(name = "Filters",
                       breaks = c("F0", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8"),
                       labels = c("F0 = Union Fusion Calls: Arriba + STAR Fusion", 
                                  "F1 = Filter Circular RNAs",
                                  "F2 = Filter by FusionAnnotator", 
                                  "F3 = Filter Read-throughs",
                                  "F4 = Collapse fusions with multiple breakpoints",
                                  "F5 = Filter by Gene Biotype: Either gene must be Protein-coding",
                                  "F6 = Create and Merge lists using various Annotations",
                                  "F7 = Expression Filter: TPM < 1 across all samples",
                                  "F8 = Tukey's Outlier Analysis: Sample is outlier for GeneA or GeneB")) +
  xlab('Filtering Steps') + ylab('Number of Fusion Calls') +
  ggtitle('Fusion Filteration Pipeline (V1)') + guides(size = F) +
  scale_y_continuous(limits = c(0, 5500), breaks = c(1000, 2000, 3000, 4000, 5000)) 
p
ggsave(filename = 'results/fusion-filtering-plot.pdf', plot = p, width = 12, height = 7)
