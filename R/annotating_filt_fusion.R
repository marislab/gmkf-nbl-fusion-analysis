# Author: Komal S Rathi
# Date: 07/19/2019
# Function: Annotate final fusions using multiple resources
# Step 3

library(reshape2)
library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)

# load fusions
load('results/fusions_v2.RData')

# A1: annotate using kinases
kinase <- read.delim('data/Kincat_Hsap.08.02.txt', stringsAsFactors = F)
total.pc.expr$kinase <- ifelse(total.pc.expr$GeneA %in% kinase$Name & total.pc.expr$GeneB %in% kinase$Name, "Both_kinase", 
                               ifelse(total.pc.expr$GeneA %in% kinase$Name, "GeneA_kinase", 
                                      ifelse(total.pc.expr$GeneB %in% kinase$Name, "GeneB_kinase", "No")))
  
# A2: annotate using TFs
tf <- read.delim('data/TRANSFAC_TFs.txt', header = F, stringsAsFactors = F)
total.pc.expr$tf <- ifelse(total.pc.expr$GeneA %in% tf$V1 & total.pc.expr$GeneB %in% tf$V1, "Both_TF", 
                           ifelse(total.pc.expr$GeneA %in% tf$V1, "GeneA_TF", 
                                  ifelse(total.pc.expr$GeneB %in% tf$V1, "GeneB_TF", "No")))

# A3: annotate oncogenes
onco <- read.delim('data/oncogenes.txt', stringsAsFactors = F)
total.pc.expr$oncogene <- ifelse(total.pc.expr$GeneA %in% onco$OncogeneName & total.pc.expr$GeneB %in% onco$OncogeneName, "Both_Onco", 
                                 ifelse(total.pc.expr$GeneA %in% onco$OncogeneName, "GeneA_Onco", 
                                        ifelse(total.pc.expr$GeneB %in% onco$OncogeneName, "GeneB_Onco", "No")))

# A4: annotate TSGs
tsg <- read.delim('data/tsgs.txt', stringsAsFactors = F)
total.pc.expr$tsg <- ifelse(total.pc.expr$GeneA %in% tsg$GeneSymbol & total.pc.expr$GeneB %in% tsg$GeneSymbol, "Both_TSG", 
                                 ifelse(total.pc.expr$GeneA %in% tsg$GeneSymbol, "GeneA_TSG", 
                                        ifelse(total.pc.expr$GeneB %in% tsg$GeneSymbol, "GeneB_TSG", "No")))



# A5: annotate TCGA fusions
tcga <- read.delim('data/tcga-fusions.txt', stringsAsFactors = F)
tcga <- separate_rows(tcga, TCGA_fusions, sep = ",")
tcga$TCGA_fusions <- gsub("_", "--", tcga$TCGA_fusions)
total.pc.expr$in_TCGA <- ifelse(total.pc.expr$Fused_Genes %in% tcga$TCGA_fusions, "TCGA_genes", "No")


# A6: annotate chimerDB
chimerDB <- read.delim('data/ChimerDB3.0_ChimerKB.txt', stringsAsFactors = F)
chimerDB <- chimerDB[,-which(colnames(chimerDB) %in% c("H_gene","H_chr","H_position","H_strand","T_gene","T_chr","T_position","T_strand"))]
colnames(chimerDB) <- paste("chimerDB", colnames(chimerDB), sep=":")
chimerDB$`chimerDB:Fusion_pair` <- sub("_","--",chimerDB$`chimerDB:Fusion_pair`)
                                       
if(length(which(chimerDB$`chimerDB:Fusion_pair` %in% total.pc.expr$Fused_Genes)) > 0){
  total.pc.expr.annotated <- merge(total.pc.expr, chimerDB, by.x = "Fused_Genes", by.y = "chimerDB:Fusion_pair", all.x = T)
} else {
  total.pc.expr.annotated <- total.pc.expr
  total.pc.expr.annotated$chimerDB <- "No"
}

# A7: PFAM domains
load('results/pfam_annotated_genes.RData')
total.pc.expr.annotated <- merge(total.pc.expr.annotated, pfam.res, by.x = 'GeneA', by.y = 'hgnc_symbol', all.x = TRUE)
colnames(total.pc.expr.annotated)[colnames(total.pc.expr.annotated) == "pfam_domain"] <- "GeneA_pfam_domain"
total.pc.expr.annotated <- merge(total.pc.expr.annotated, pfam.res, by.x = 'GeneB', by.y = 'hgnc_symbol', all.x = TRUE)
colnames(total.pc.expr.annotated)[colnames(total.pc.expr.annotated) == "pfam_domain"] <- "GeneB_pfam_domain"
total.pc.expr.annotated[is.na(total.pc.expr.annotated)] <- "No"

# save fusions
save(total.pc.expr.annotated, file = 'results/fusions_v3.RData')
