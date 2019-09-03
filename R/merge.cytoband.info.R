merge.cytoband.info <- function(fusion.calls, cytoband){
  
  # format cytoband info
  cytoband$V4 <- gsub('[0-9.]','',cytoband$V4)
  cytoband <- cytoband[grep("_", cytoband$V1, invert = TRUE),]
  cytoband$V1 <- gsub("chr", "", cytoband$V1)
  cytoband <- cytoband %>% group_by(V1, V4) %>%
    filter(V1 != "M") %>%
    summarise(V2 = min(V2),
              V3 = max(V3)) %>%
    as.data.frame()
  
  # add cytoband info to chromosome
  tmp <- unique(fusion.calls[,c('GeneA_bp','GeneB_bp')])
  tmp <- cbind(tmp, colsplit(tmp$GeneA_bp,':', names = c("chrA","startA")))
  tmp$endA <- tmp$startA
  tmp <- cbind(tmp, colsplit(tmp$GeneB_bp,':', names = c("chrB","startB")))
  tmp$endB <- tmp$startB
  one <- unique(tmp[,c('chrA','startA','endA','GeneA_bp')])
  two <- unique(tmp[,c('chrB','startB','endB','GeneB_bp')])
  
  cytoband <- cytoband[,c(1,3,4,2)]
  cytoband$V4 <- paste0(cytoband$V1,cytoband$V4)
  colnames(cytoband) <- c("chr","start","end","name")
  
  # for chrA
  subject <- with(cytoband, GRanges(chr, IRanges(start = start, end = end, names = name)))
  query <- with(one, GRanges(chrA, IRanges(start = startA, end = endA, names = GeneA_bp)))
  res <- findOverlaps(query = query, subject = subject, type = "within")
  res.df <- data.frame(one[queryHits(res),], cytoband[subjectHits(res),])
  res.geneA <- res.df[,c("GeneA_bp","name")]
  colnames(res.geneA)[2] <- "GeneA_cyto"
  
  # for chrB
  query <- with(two, GRanges(chrB, IRanges(start = startB, end = endB, names = GeneB_bp)))
  res <- findOverlaps(query = query, subject = subject, type = "within")
  res.df <- data.frame(two[queryHits(res),], cytoband[subjectHits(res),])
  res.geneB <- res.df[,c("GeneB_bp","name")]
  colnames(res.geneB)[2] <- "GeneB_cyto"
  
  # merge with fusion calls
  fusion.calls <- merge(fusion.calls, res.geneA, by = 'GeneA_bp')
  fusion.calls <- merge(fusion.calls, res.geneB, by = 'GeneB_bp')
  nums <- c(1:22, "X", "Y")
  cyto <- c("p", "q")
  levels <- paste(rep(nums, each = length(cyto)), cyto, sep = "")
  fusion.calls$GeneA_cyto <- factor(fusion.calls$GeneA_cyto, levels = levels)
  fusion.calls$GeneB_cyto <- factor(fusion.calls$GeneB_cyto, levels = levels)
  
  return(fusion.calls)
}