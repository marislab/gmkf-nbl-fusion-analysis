# format filtered Arriba file to run fusion annotator

setwd('~/Projects/Maris-lab/GMKF_NBL/')
library(dplyr)
library(tidyr)

source('R/format.arriba.R')
load('data/arriba_merged.RData')
ar <- format.arriba(ar = arriba.merge)
head(ar)
fus <- unique(ar$FusionName)
write.table(fus, file = 'results/arriba_fusion_fusionannot_input.txt', quote = F, row.names = F, col.names = F)

# # format annotation
# load('data/expr.RData')
# expr.mat <- expr.mat[grep("^PAR_Y", expr.mat$gene_symbol, invert = T),]
# expr.mat$tpm_means <- apply(expr.mat[,3:ncol(expr.mat)], MARGIN = 1, FUN = mean)
# expr.mat <- expr.mat[order(expr.mat$gene_symbol, expr.mat$tpm_means, decreasing = TRUE),]
# expr.mat <- expr.mat[!duplicated(expr.mat$gene_symbol),]
# gene.ids <- expr.mat$gene_id
# rm(expr.mat)
# 
# annot <- read.delim('data/gencode.v27.primary_assembly.annotation.txt', stringsAsFactors = F)
# annot <- annot[grep("PAR", annot$gene_id, invert = TRUE),]
# annot <- annot[which(annot$gene_id %in% gene.ids),]
# annot$gene_id <- gsub("[.].*", "", annot$gene_id)

# # merge arriba with annotation to get ensembl ids
# ar <- merge(ar, annot, by.x = 'gene1', by.y = 'gene_symbol', all.x = TRUE)
# colnames(ar)[colnames(ar) == "gene_id"] <- "5p_ensembl"
# ar <- merge(ar, annot, by.x = 'gene2', by.y = 'gene_symbol', all.x = TRUE)
# colnames(ar)[colnames(ar) == "gene_id"] <- "3p_ensembl"


