# Author: Komal S Rathi
# Date: 08/15/2019
# Function to format expression data for Diff expression analysis with GTEx tissues
# compared to Adult tumors, a list of genes that are differentially expressed in NBL tumors

# create a matrix with gene ids are rownames
# meta data with 'study_id','disease','sample_id'
load('data/expr_counts.RData')
expr.mat <- expr.mat[grep("^PAR_Y", expr.mat$gene_symbol, invert = TRUE),]
expr.mat$gene_id <- gsub("[.].*", "", expr.mat$gene_id)
rownames(expr.mat) <- expr.mat$gene_id
expr.mat <- expr.mat[,3:ncol(expr.mat)]
tumor.mat <- expr.mat
tumor.mat <- tumor.mat[,order(colnames(tumor.mat))]

load('data/clin_mapped.RData')
clin$study_id <- "GMKF_NBL"
clin$sample_id <- clin$name
clin$disease <- "GMKF_NBL"
clin <- clin[,c('study_id','disease','sample_id')]
rownames(clin) <- clin$sample_id
clin <- clin[order(clin$sample_id),]

if(identical(colnames(tumor.mat), rownames(clin))){
  print("identical")
  tumormeta <- clin
}  

# add annotation
annotation <- read.delim('data/gencode.v27.primary_assembly.annotation.txt')
annotation$biotype <- NULL
annotation <- annotation[grep('_PAR_Y$', annotation$gene_id, invert = T),]
annotation$gene_id <- gsub('[.].*','',annotation$gene_id)

save(annotation, tumor.mat, tumormeta, file = 'results/expr_counts_forDE.RData')
