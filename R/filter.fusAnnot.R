# function to filter red herrings

filter.fusAnnot <- function(caller, red.herring, labels.to.keep){
  
  # FusionAnnotator filters
  caller.FusAnnot <- unique(caller[grep(red.herrings, caller$annots, ignore.case = TRUE),])
  caller.FusAnnot.rm <- unique(caller.FusAnnot[grep(labels.to.keep, caller.FusAnnot$annots, ignore.case = TRUE),])
  caller.FusAnnot.rm <- unique(caller.FusAnnot.rm[which(caller.FusAnnot.rm$Confidence == "high"),"Fused_Genes"])
  caller.FusAnnot <- caller.FusAnnot[-which(caller.FusAnnot$Fused_Genes %in% caller.FusAnnot.rm),]
  caller.FusAnnot <- unique(caller.FusAnnot$Fused_Genes)
  caller.FusAnnot.rev <- unique(unlist(lapply(strsplit(caller.FusAnnot, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
  fusions.to.remove <- unique(c(caller.FusAnnot, caller.FusAnnot.rev))
  
  # # arriba FusionAnnotator
  # arriba.FusAnnot <- unique(ar[grep(red.herrings, ar$annots, ignore.case = TRUE),])
  # arriba.FusAnnot.rm <- unique(arriba.FusAnnot[grep(labels.to.keep, arriba.FusAnnot$annots, ignore.case = TRUE),])
  # arriba.FusAnnot.rm <- unique(arriba.FusAnnot.rm[which(arriba.FusAnnot.rm$Confidence == "high"),"FusionName"])
  # arriba.FusAnnot <- arriba.FusAnnot[-which(arriba.FusAnnot$FusionName %in% arriba.FusAnnot.rm),]
  # arriba.FusAnnot <- unique(arriba.FusAnnot$FusionName)
  # 
  # # StarFusion FusionAnnotator
  # sf.FusAnnot <- unique(sf[grep(red.herrings, sf$annots, ignore.case = TRUE),])
  # sf.FusAnnot.rm <- unique(sf.FusAnnot[grep(labels.to.keep, sf.FusAnnot$annots, ignore.case = TRUE),])
  # sf.FusAnnot.rm <- unique(sf.FusAnnot.rm[which(sf.FusAnnot.rm$Confidence == "high"),"FusionName"])
  # sf.FusAnnot <- sf.FusAnnot[-which(sf.FusAnnot$FusionName %in% sf.FusAnnot.rm),]
  # sf.FusAnnot <- unique(sf.FusAnnot$FusionName)
  # 
  # fusions.to.remove <- unique(c(sf.FusAnnot, arriba.FusAnnot))
  # fusions.to.remove.rev <- unique(unlist(lapply(strsplit(fusions.to.remove, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
  # fusions.to.remove <- unique(c(fusions.to.remove, fusions.to.remove.rev))
  return(fusions.to.remove)
}
