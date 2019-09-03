# Author: Komal S. Rathi
# Date: 08/30/2019
# Function to format and filter Manta RNA-seq Str. Variants

setwd('~/Projects/Maris-lab/GMKF_NBL/')
load('data/str_vars_merged.RData')

str.var.mat <- str.var.mat[grep('chr', str.var.mat$CHROM),]
str.var.mat <- str.var.mat[grep('chr', str.var.mat$ALT),]
str.var.mat$ALT <- gsub(".*chr","",str.var.mat$ALT)
str.var.mat$ALT <- gsub("[[]|[]]|[ATGC]","",str.var.mat$ALT)
str.var.mat$CHROM <- gsub("chr","", str.var.mat$CHROM)
str.var.mat$GeneA_bp <- paste0(str.var.mat$CHROM,':',str.var.mat$POS)
colnames(str.var.mat)[colnames(str.var.mat) == "ALT"] <- "GeneB_bp"
str.var.mat <- unique(str.var.mat[,c("GeneA_bp","GeneB_bp","Sample")])
save(str.var.mat, file = 'data/str_vars_merged_short.RData')
