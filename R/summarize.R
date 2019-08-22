# Author: Komal S Rathi
# Date: 07/19/2019
# Function: Summarize and visualize final fusions 

load('results/fusions_v3.RData')

# add differential expression results to fusion data
# gtex tissue comparison
degs <- readRDS('results/GMKF_NBL_PST_224.RDS')
degs <- degs[degs$maxAdjPVal < 0.05,]
degs <- degs[abs(degs$minLFC) > 1,] # 7019 differentially expressed
degs$deg_annot[degs$minLFC > 1 & degs$maxAdjPVal < 0.05] <- "Up"
degs$deg_annot[degs$minLFC < -1 & degs$maxAdjPVal < 0.05] <- "Down"
degs <- degs[,c('gene_id','gene_symbol','maxAdjPVal','minLFC','deg_annot')]

# remove duplicate gene ids
ids <- read.delim('results/geneids_to_keep.txt')
ids <- gsub('[.].*','',ids$x)
degs <- degs[which(degs$gene_id %in% ids),]

# annotate
total.pc.expr.annotated <- merge(total.pc.expr.annotated, degs[,c('gene_symbol','deg_annot')], by.x = 'GeneA', by.y = 'gene_symbol', all.x = TRUE)
colnames(total.pc.expr.annotated)[colnames(total.pc.expr.annotated) == "deg_annot"] <- "GeneA_expr_status"
total.pc.expr.annotated <- merge(total.pc.expr.annotated, degs[,c('gene_symbol','deg_annot')], by.x = 'GeneB', by.y = 'gene_symbol', all.x = TRUE)
colnames(total.pc.expr.annotated)[colnames(total.pc.expr.annotated) == "deg_annot"] <- "GeneB_expr_status"

# remove all fusions where expression does not change compared to normal dataset
tmp <- total.pc.expr.annotated[!is.na(total.pc.expr.annotated$GeneA_expr_status) | 
                                 !is.na(total.pc.expr.annotated$GeneB_expr_status),]
