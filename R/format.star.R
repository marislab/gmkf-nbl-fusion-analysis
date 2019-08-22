# Author: Komal S. Rathi
# Date: 08/13/2019
# Function: Function to format and filter star output

format.star <- function(sf){
  sf <- sf[-which(sf$SpanningFragCount-sf$JunctionReadCount > 10 |
                    sf$JunctionReadCount == 0 |
                    sf$LargeAnchorSupport == "NO_LDAS"),]
  sf$LeftBreakpoint <- gsub('^chr', '', sf$LeftBreakpoint)
  sf$RightBreakpoint <- gsub('^chr', '', sf$RightBreakpoint)
  colnames(sf)[c(6,8)] <- c('Gene1_pos','Gene2_pos')
  sf$Fusion_Type <- "NA"
  sf$Fusion_Type[which(sf$PROT_FUSION_TYPE == "INFRAME")] <- 'in-frame'
  sf$Fusion_Type[grep("FRAMESHIFT",sf$PROT_FUSION_TYPE)] <- 'frameshift'
  sf$Fusion_Type[-which(sf$PROT_FUSION_TYPE %in% c("INFRAME","FRAMESHIFT"))] <- 'other'
  sf$Caller <- 'STARFusion'
  colnames(sf)[1] <- 'FusionName'
  sf$Confidence <- "NA"
  colnames(sf)[colnames(sf) == "sample_name"] <- "Sample"
  return(sf)
}