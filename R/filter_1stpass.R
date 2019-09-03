# script to filter low confidence calls from Arriba

setwd('~/Projects/Maris-lab/GMKF_NBL/')

# low confidence calls
load('data/arriba_merged.RData')
# arriba.merge <- arriba.merge[which(arriba.merge$confidence %in% c("low", "medium")),]
arriba.merge <- unique(arriba.merge[,c("breakpoint1","breakpoint2","sample_name","confidence")])

# str variants
load('data/str_vars_merged_short.RData')
res <- arriba.merge[(arriba.merge$breakpoint1 %in% str.var.mat$GeneA_bp & 
                       arriba.merge$sample_name %in% str.var.mat$Sample) |
               (arriba.merge$breakpoint2 %in% str.var.mat$GeneB_bp & 
                  arriba.merge$sample_name %in% str.var.mat$Sample),]

