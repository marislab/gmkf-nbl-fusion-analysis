# Author: Komal S. Rathi
# Date: 08/19/2019
# Function: Create Summary Plots
source('~/Projects/DGD_Mendelian_RNASeq/R/pubTheme.R')

# Line plot
df <- data.frame(Filter = c("F0","F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8"), 
                 Label = c("Combined Calls", "Filter GTEx Recurrent", "Filter Read-throughs", 
                           "Filter Circular RNAs", "Filter Duplicates", "Filter by Gene Biotype", 
                           "Annotation filters", "Expression Filter", "Outlier Analysis"),
                 Val = c(4926, 4926, 1886, 1850, 1650, 1461, 520, 520, 171))
df$Sample <- "GMKF_NBL"


p <- ggplot(df, aes(x = Filter, y = Val, group = Sample)) + 
  geom_point(aes(color = Filter), pch = 18, size = 6) + 
  geom_line(linetype = "dotted") + 
  geom_text(aes(label = Val, hjust = 0.5, vjust = -1)) + theme_Publication3() + 
  theme(legend.position = "right") +
  scale_color_discrete(name = "Filters",
                       breaks = c("F0", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8"),
                       labels = c("F0 = Union Fusion Calls: Arriba + STAR Fusion", 
                                  "F1 = Filter GTEx Recurrent",
                                  "F2 = Filter Read-throughs", 
                                  "F3 = Filter Circular RNAs",
                                  "F4 = Reduce duplicated fusions with multiple breakpoints",
                                  "F5 = Filter by Gene Biotype: Either gene must be Protein-coding",
                                  "F6 = Create and Merge lists",
                                  "F7 = Expression Filter: TPM < 1 across all samples",
                                  "F8 = Outlier Analysis: Sample is outlier for GeneA or GeneB")) +
  xlab('Filtering Steps') + ylab('Number of Fusion Calls') +
  ggtitle('Fusion Filteration Pipeline (V1)') + guides(size = F) +
  scale_y_continuous(limits = c(0, 5500), breaks = c(1000, 2000, 3000, 4000, 5000)) 
p
ggsave(filename = 'results/fusion-filtering-plot.pdf', plot = p, width = 12, height = 7)
