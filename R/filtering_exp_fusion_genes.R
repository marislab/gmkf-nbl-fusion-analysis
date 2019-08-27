# Author: Komal S Rathi
# Date: 07/19/2019
# Function: Filter fusion calls based on expression
# Remove fusion calls, where both fusion genes are not expressed
# Step 2

library(reshape2)
library(data.table)
library(plyr)

source('R/detect.outliers.R')

# Load and Format expression data matrix
# collapse to gene symbols
load('data/expr.RData')
expr.mat <- expr.mat[grep('^PAR_Y_', expr.mat$gene_symbol, invert = TRUE),]
expr.mat$tpm_means <- apply(expr.mat[,3:ncol(expr.mat)], MARGIN = 1, FUN = mean)
expr.mat <- expr.mat[order(expr.mat$gene_symbol, expr.mat$tpm_means, decreasing = TRUE),]
expr.mat <- expr.mat[!duplicated(expr.mat$gene_symbol),]
write.table(expr.mat$gene_id, file = 'results/geneids_to_keep.txt', quote = F, sep = "\t", row.names = F)

# add filter for < 1 FPKM
expr.mat$expressed <- !(apply(expr.mat[,3:226], 1, FUN = function(x) all(x < 1)))
for.quantiles <- expr.mat
expr.mat <- expr.mat[,c('gene_symbol','tpm_means','expressed')]

# annotate to fusion file
# only keep fusions where both genes are expressed
# (n = 628)
load('results/fusions_v1.RData')
total <- cbind(total, colsplit(total$Fused_Genes, '--', names = c("GeneA", "GeneB")))
total <- merge(total, expr.mat, by.x = 'GeneA', by.y = 'gene_symbol')
colnames(total)[colnames(total) %in% c("tpm_means","expressed")] <- paste0('GeneA_', c("tpm_means","expressed"))
total <- merge(total, expr.mat, by.x = 'GeneB', by.y = 'gene_symbol')
colnames(total)[colnames(total) %in% c("tpm_means","expressed")] <- paste0('GeneB_', c("tpm_means","expressed"))
total.expr <- total[which(total$GeneA_expressed == TRUE | total$GeneB_expressed == TRUE),]

# calculate quantiles and add to the fusion table as a 1-0 binary
if(file.exists('data/expr_quantile_change.RData')){
  print("exists")
  load('data/expr_quantile_change.RData')
} else {
  for.quantiles$expressed <- NULL
  for.quantiles$tpm_means <- NULL
  rownames(for.quantiles) <- for.quantiles$gene_symbol
  for.quantiles.new <- apply(for.quantiles[,3:ncol(for.quantiles)], 1, detect.outliers)
  for.quantiles.new <- as.data.frame(t(for.quantiles.new))
  for.quantiles.new <- melt(as.matrix(for.quantiles.new))
  colnames(for.quantiles.new) <- c("gene_symbol","sample_id","expression_change")
  expr.quantiles <- for.quantiles.new
  
  # save to file
  save(expr.quantiles, file = 'data/expr_quantile_change.RData')
}

# annotate fusions with relative expression change (Expression outliers)
# (n = 171)
total.pc.expr <- total.expr
total.pc.expr <- merge(total.pc.expr, expr.quantiles, by.x = c("Sample","GeneA"), by.y = c("sample_id","gene_symbol"), all.x = TRUE)
colnames(total.pc.expr)[colnames(total.pc.expr) == "expression_change"] <- "GeneA_expression_change"
total.pc.expr <- merge(total.pc.expr, expr.quantiles, by.x = c("Sample","GeneB"), by.y = c("sample_id","gene_symbol"), all.x = TRUE)
colnames(total.pc.expr)[colnames(total.pc.expr) == "expression_change"] <- "GeneB_expression_change"
total.pc.expr <- total.pc.expr[which(total.pc.expr$GeneA_expression_change != 0 | total.pc.expr$GeneB_expression_change != 0),]

# now redefine note column based on the expression filters applied
# need to only define the recurrent fusions and recurrently fused genes
total.pc.expr$note <- NULL
total.pc.expr <- unique(total.pc.expr)
l1 <- read.delim('results/inframe_fusions_L1.txt')
l2 <- read.delim('results/Onco_TSG_TCGA_fusions_L2.txt')
l3 <- read.delim('results/BothCallers_fusions_L3.txt')
l4 <- read.delim('results/Recurrent_Fusions_L4.txt')
l5 <- read.delim('results/Recurrently_Fused_Genes_L5.txt')
total.pc.expr$In_frame <- ifelse(total.pc.expr$Sample %in% l1$Sample &
                                   total.pc.expr$Fused_Genes %in% l1$Fused_Genes &
                                   total.pc.expr$Caller %in% l1$Caller &
                                   total.pc.expr$Fusion_Type %in% l1$Fusion_Type, "Yes", "No")
total.pc.expr$onco_tsg_tcga <- ifelse(total.pc.expr$Sample %in% l2$Sample &
                                        total.pc.expr$Fused_Genes %in% l2$Fused_Genes &
                                        total.pc.expr$Caller %in% l2$Caller &
                                        total.pc.expr$Fusion_Type %in% l2$Fusion_Type, "Yes", "No")
total.pc.expr$both_callers <- ifelse(total.pc.expr$Sample %in% l3$Sample &
                                       total.pc.expr$Fused_Genes %in% l3$Fused_Genes &
                                       total.pc.expr$Caller %in% l3$Caller &
                                       total.pc.expr$Fusion_Type %in% l3$Fusion_Type, "Yes", "No")
total.pc.expr$rec_fusions <- ifelse(total.pc.expr$Sample %in% l4$Sample &
                                      total.pc.expr$Fused_Genes %in% l4$Fused_Genes &
                                      total.pc.expr$Caller %in% l4$Caller &
                                      total.pc.expr$Fusion_Type %in% l4$Fusion_Type, "Yes", "No")
total.pc.expr$rec_fused_genes <- ifelse(total.pc.expr$Sample %in% l5$Sample &
                                          total.pc.expr$Fused_Genes %in% l5$Fused_Genes &
                                          total.pc.expr$Caller %in% l5$Caller &
                                          total.pc.expr$Fusion_Type %in% l5$Fusion_Type, "Yes", "No")
plyr::count(total.pc.expr$In_frame)
plyr::count(total.pc.expr$onco_tsg_tcga)
plyr::count(total.pc.expr$both_callers)
plyr::count(total.pc.expr$rec_fusions)
plyr::count(total.pc.expr$rec_fused_genes)

# save
save(total.pc.expr, file =  'results/fusions_v2.RData')

