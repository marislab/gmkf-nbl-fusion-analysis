# Author: Komal S. Rathi
# Date: 08/13/2019
# Function: Add biotypes to gene symbols 
# Filter fusions where none of the genes is protein_coding

add.biotypes <- function(df, annot){
  annot <- annot[,c("gene_symbol", "biotype")]
  df <- merge(df, annot, by.x = "GeneA", by.y = "gene_symbol")
  colnames(df)[colnames(df) == "biotype"] <- "GeneA.biotype"
  df <- merge(df, annot, by.x = "GeneB", by.y = "gene_symbol")
  colnames(df)[colnames(df) == "biotype"] <- "GeneB.biotype"
  
  # Filter F1: either gene should be protein coding
  df <- df[which(df$GeneA.biotype == "protein_coding" | df$GeneB.biotype == "protein_coding"),]
  return(df)
}
