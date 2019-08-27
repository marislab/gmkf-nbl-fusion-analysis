# Author: Komal S. Rathi
# Date: 08/13/2019
# Function to format and filter arriba output
# Add output of fusion annotator to arriba output

format.arriba <- function(ar, arriba.FusAnnot){
  
  colnames(ar)[1] <- "gene1"
  
  # # F1: filter low confidence intergenic fusions
  # ar <- ar[-which((ar$site1 == "intergenic" | ar$site2 == "intergenic") & 
  #                   (ar$confidence == "low")),]
  
  # expand intergenic fusions into new rows
  ar <- ar %>%
    mutate(gene1 = strsplit(as.character(gene1), ",")) %>%
    unnest(gene1) %>%
    as.data.frame()
  ar <- ar %>%
    mutate(gene2 = strsplit(as.character(gene2), ",")) %>%
    unnest(gene2) %>%
    as.data.frame()
  ar$gene1 <- gsub('[(].*','',ar$gene1)
  ar$gene2 <- gsub('[(].*','',ar$gene2)
  
  # F2: filter out Circular RNAs
  ar <- ar[grep('duplication/non-canonical_splicing', ar$type, invert = TRUE),]
  
  # F3: filter truncated genes/fusions
  # genes fused head-on or tail-tail cannot yeild a chimeric protein
  # since one of the genes is transcribed from the wrong strand
  ar <- ar[grep("5\'-5\'|3\'-3\'", ar$type, invert = TRUE),]
  
  # F4: filter false positives through supporting reads
  ar <- ar[-(which(ar$discordant_mates-(ar$split_reads1+ar$split_reads2) > 10 |
                     ar$split_reads1 + ar$split_reads2 == 0)),]
  
  # format
  ar$LeftBreakpoint <- gsub('^chr','',ar$breakpoint1)
  ar$RightBreakpoint <- gsub('^chr','',ar$breakpoint2)
  ar$Fusion_Type <- ar$reading_frame
  ar$Fusion_Type[grep("out-of-frame",ar$Fusion_Type)] <- "frameshift"
  ar$Fusion_Type[-which(ar$Fusion_Type %in% c("in-frame","frameshift"))] <- 'other'
  ar$Caller <- 'arriba'
  colnames(ar)[colnames(ar) == "sample_name"] <- "Sample"
  ar$FusionName <- paste0(gsub(",","/",ar$gene1),"--",gsub(",","/",ar$gene2))
  ar$SpanningFragCount <- ar$discordant_mates
  ar$JunctionReadCount <- ar$split_reads1+ar$split_reads2
  ar$Confidence <- ar$confidence
  
  ar$GeneA_bp <- ar$breakpoint1
  ar$GeneB_bp <- ar$breakpoint2
  
  # add fusion annotator output
  ar <- merge(ar, arriba.FusAnnot, by = 'FusionName')
  
  # add gene names
  ar <- cbind(ar, colsplit(ar$FusionName, '--', c("GeneA", "GeneB")))

  return(ar)
}
