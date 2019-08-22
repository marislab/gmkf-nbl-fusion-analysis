# Author: Komal S Rathi
# Date: 07/19/2019
# Function: Filter fusion calls
# Step 1

library(reshape2)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
detach("package:plyr", unload=TRUE)

# source functions
source('R/format.arriba.R')
source('R/format.star.R')
source('R/add.biotypes.R')

# Step 1
# load data
load('data/arriba_merged.RData')
load('data/star_merged.RData')
sf <- as.data.frame(star.merge)
ar <- as.data.frame(arriba.merge)
rm(star.merge, arriba.merge)

# star fusion
sf <- format.star(sf)
sf.total <- unique(sf[,c('FusionName','Sample','Caller','Fusion_Type','JunctionReadCount','SpanningFragCount','Confidence')])

# arriba fusion
# format arriba output for intergenic fusions
ar <- format.arriba(ar)
ar.total <- unique(ar[,c('FusionName','Sample','Caller','Fusion_Type','JunctionReadCount','SpanningFragCount','Confidence')])

# identify readthrough fusions
arriba.FusAnnot <- read.delim('results/arriba_fusions_annotated.txt', header = F, stringsAsFactors = F)
arriba.FusAnnot <- unique(arriba.FusAnnot[grep("readthrough|NEIGHBORS|GTEx_Recurrent", arriba.FusAnnot$V2),'V1'])
sf.rt <- unique(sf[grep('readthrough|NEIGHBORS|GTEx_Recurrent', sf$annots), 'FusionName'])
ar.rt <- unique(ar[grep("read-through|non-canonical_splicing", ar$type), 'FusionName'])
rts <- unique(c(sf.rt, ar.rt, arriba.FusAnnot))
rts.rev <- unique(unlist(lapply(strsplit(rts, '--'), FUN = function(x) paste0(x[2],'--',x[1]))))
rts <- unique(c(rts, rts.rev))
rm(sf.rt, ar.rt, arriba.FusAnnot)

# merge both callers
all.callers <- rbind(ar.total, sf.total)
rm(ar, sf, ar.total, sf.total)

# histology (no sub-histology for NBL)
load('data/clin_mapped.RData')

# merge callers and clinical information
all.callers <- merge(all.callers, clin, by.x = "Sample", by.y = "name")
all.callers <- unique(all.callers)
colnames(all.callers)[2] <- "Fused_Genes"

# Filter F0: remove all non-expressed genes
# all.callers <- filter.exp(all.callers)

# Filter F1: remove read-throughs (n = 1886)
all.callers <- all.callers[-which(all.callers$Fused_Genes %in% rts),]

# Filter F2: filter Circular RNAs (n = 1850)
all.callers <- cbind(all.callers, colsplit(all.callers$Fused_Genes, pattern = '--', names = c("GeneA","GeneB")))
# all.callers <- all.callers[!(all.callers$GeneA == all.callers$GeneB & 
#                         all.callers$Confidence != "high" & 
#                         all.callers$Fusion_Type != "in-frame"),]
circ.rna <- all.callers[which(all.callers$GeneA == all.callers$GeneB),]
write.table(circ.rna, file = 'results/circ_RNA_txt', quote = F, sep = "\t", row.names = F)
all.callers <- all.callers[which(all.callers$GeneA != all.callers$GeneB),]

# Filter F3: reduce duplications (n = 1650)
# reduce duplicate fusions to keep one fusion per caller with the highest JunctionReadCount and Confidence
# assign 3 to high confidence, 2 is medium and 1 is low
all.callers$Confidence <- ifelse(all.callers$Confidence %in% c("high", "NA"), 3, ifelse(all.callers$Confidence == "medium", 2, 1))
all.callers <- all.callers %>%
  group_by_at(vars(-JunctionReadCount, -SpanningFragCount, -Confidence)) %>%
  mutate(max.junc = ifelse(JunctionReadCount == max(JunctionReadCount), 'yes', 'no'),
         max.frag = ifelse(SpanningFragCount == max(SpanningFragCount), 'yes', 'no'),
         max.conf = ifelse(Confidence == max(Confidence), 'yes', 'no')) %>%
  filter(max.junc == "yes" & max.frag == "yes" & max.conf == "yes") %>%
  as.data.frame()

# Filter F4: annotate gene1 and gene2 with biotype
# filter any fusion where none of the genes is protein coding
# (n = 1461)
# load('results/fusions_v1.RData')
annot <- read.delim('data/gencode.v27.primary_assembly.annotation.txt', stringsAsFactors = F)
annot <- annot %>% group_by(gene_symbol) %>%
  summarise(gene_id = toString(unique(gene_id)),
            biotype = toString(unique(biotype))) %>%
  unique() %>%
  as.data.frame()
all.callers <- add.biotypes(df = all.callers, annot = annot)

# List L1: keep all in-frame fusions
# (n = 125)
inframe.fus <- all.callers %>% 
  filter(Fusion_Type == "in-frame") %>%
  as.data.frame()
inframe.fus$note <- "In-frame fusions"
tmp <- unique(inframe.fus[,c('Fused_Genes','Sample','Caller','Fusion_Type')])
write.table(tmp, file = 'results/inframe_fusions_L1.txt', quote = F, sep = "\t", row.names = F)
  
# List L2: list of fusions to be added irrespective of frame
# driver fusions from Literature Jo Lynne
# no driver fusions for Neuroblastoma
# lit.genes <- read.delim('../../references/driver-fusions-v4.txt', stringsAsFactors = F)
# lit.genes <- lit.genes[!is.na(lit.genes$FusionPartner),]
# modify input gene lists
tsgs <- read.delim('data/tsgs.txt', stringsAsFactors = F)
tsgs <- unique(tsgs$GeneSymbol)
onco <- read.delim('data/oncogenes.txt', stringsAsFactors = F)
onco <- unique(onco$OncogeneName)
tcga <- read.delim('data/tcga-fusions.txt', stringsAsFactors = F)
tcga <- unique(unlist(strsplit(tcga$TCGA_fusions,",")))
fusions.to.search <- gsub("_","--",tcga)

# only run if you have a list of lit genes
if(exists('lit.genes')) {
  for(i in 1:nrow(lit.genes)){
    genes.to.search <- lit.genes[i,2]
    fusions.to.search<-gsub("_","--",lit.genes[i,3])
    genes.to.search<-unlist(strsplit(genes.to.search, ','))
    print(fusions.to.search)
    genes.to.search <- c(paste0('^',genes.to.search,'-'), paste0('-',genes.to.search,'$'))
    print(genes.to.search)
    if(fusions.to.search == ""){
      print("no fusions to check")
    } else {
      fusions.to.search <- paste0('^',fusions.to.search,'$')
      genes.to.search <- c(genes.to.search, fusions.to.search)
    }
    genes.to.search <- paste0(genes.to.search, collapse = '|')
    hist.to.search <- lit.genes[1,1]
    getfusions <- all.callers[grep(genes.to.search, all.callers$Fused_Genes),]
    getfusions <- getfusions[which(getfusions$Histology.Broad %in% hist.to.search),]
    getfusions <- unique(getfusions)
    if(nrow(getfusions) == 0){
      print(hist.to.search)
      print(genes.to.search)
    }
    if(i == 1){
      lit.genes.add <- getfusions
    } else {
      lit.genes.add <- rbind(lit.genes.add, getfusions)
    }
  }
  lit.genes.add <- unique(lit.genes.add)
  colnames(lit.genes.add) <- colnames(all.callers)
  print("Literature gene fusions")
  lit.genes.add$note <- "known oncogenic fusion"
}

# tsgs,onco,tcga fusion genes 
# (n = 299)
genes.to.search <- unique(c(tsgs, onco))
genes.to.search <- c(paste0('^',genes.to.search,'--'), paste0('--',genes.to.search,'$'))
tsgs_onco <- all.callers[unlist(lapply(genes.to.search, function(x) grep(x, all.callers$Fused_Genes))),]
fusion <- all.callers[unlist(lapply(fusions.to.search, function(x) grep(x, all.callers$Fused_Genes))),]
tsgs_onco_fusion <- rbind(tsgs_onco, fusion)
tsgs_onco_fusion$note <- "Found in onco/tsgs/tcga list"
if(exists('lit.genes.add')){
  to.add <- rbind(lit.genes.add, tsgs_onco_fusion)
} else {
  to.add <- tsgs_onco_fusion
}
tmp <- unique(to.add[,c('Fused_Genes','Sample','Caller','Fusion_Type')])
write.table(tmp, file = 'results/Onco_TSG_TCGA_fusions_L2.txt', quote = F, sep = "\t", row.names = F)


# List L3: Found by two callers and Fusion Type is well annotated
# (n = 52)
two.callers <- all.callers %>%
  as.data.frame() %>%
  filter(Fusion_Type != "other") %>%
  group_by(Sample, Fused_Genes, Fusion_Type) %>%
  unique() %>%
  mutate(caller = toString(Caller), caller.count = n()) %>%
  filter(caller.count >= 2) %>%
  unique() %>%
  as.data.frame()
two.callers$note <- "All callers"
tmp <- unique(two.callers[,c('Fused_Genes','Sample','Caller','Fusion_Type')])
write.table(tmp, file = 'results/BothCallers_fusions_L3.txt', quote = F, sep = "\t", row.names = F)

# List L4: Found in at least 2 samples of the same histology and Fusion Type is well annotated
# For recurrent fusions, filter out low confidence calls
# (n = 32)
rec.fusions <- all.callers %>%
  as.data.frame() %>%
  filter(Fusion_Type != "other" & Confidence > 1) %>%
  group_by(Fused_Genes, Caller, Fusion_Type) %>%
  mutate(sample.count = n()) %>%
  filter(sample.count > 1) %>%
  unique() %>%
  as.data.frame()
rec.fusions$note <- "Recurrent fusions"
tmp <- unique(rec.fusions[,c('Fused_Genes','Sample','Caller','Fusion_Type')])
write.table(tmp, file = "results/Recurrent_Fusions_L4.txt", quote = F, sep = "\t", row.names = F)

# List L5: GeneB or GeneA gene recurrently fused within a histology (>= 5 genes)
# For recurrently fused genes, filter out low confidence calls
# (n = 130)
# rec <- cbind(all.callers, colsplit(all.callers$Fused_Genes, pattern = '--', names = c("GeneA","GeneB")))
all.callers$Histology.Broad <- "Neuroblastoma"
rec <- all.callers
# (n = 32)
rec.geneA <- rec %>%
  filter(Confidence > 1) %>%
  group_by(Histology.Broad, Fused_Genes, GeneA, Caller, Fusion_Type) %>%
  dplyr::select(Histology.Broad, Fused_Genes, GeneA, Caller, Fusion_Type) %>%
  unique() %>% 
  group_by(Histology.Broad, GeneA, Caller, Fusion_Type) %>%
  summarise(GeneA.ct = n()) %>%
  filter(GeneA.ct >= 3) %>% as.data.frame()

# (n = 8)
rec.geneB <- rec %>%
  filter(Confidence > 1) %>%
  group_by(Histology.Broad, Fused_Genes, GeneB, Caller, Fusion_Type) %>%
  dplyr::select(Histology.Broad, Fused_Genes, GeneB, Caller, Fusion_Type) %>%
  unique() %>%
  group_by(Histology.Broad, GeneB, Caller, Fusion_Type) %>%
  summarise(GeneB.ct = n()) %>%
  filter(GeneB.ct >= 3) %>% as.data.frame()
  
rec.geneA <- merge(rec.geneA, rec)
rec.geneA$note <- "GeneA recurrently fused"
rec.geneB <- merge(rec.geneB, rec)
rec.geneB$note <- "GeneB recurrently fused"
rec.geneA <- unique(rec.geneA[,colnames(to.add)])
rec.geneB <- unique(rec.geneB[,colnames(to.add)])
rec.fused.genes <- unique(rbind(rec.geneA, rec.geneB))
tmp <- unique(rec.fused.genes[,c('Fused_Genes','Sample','Caller','Fusion_Type')])
write.table(tmp, file = "results/Recurrently_Fused_Genes_L5.txt", quote = F, sep = "\t", row.names = F)

# merge all five lists (n = 520)
# collapse by note
cols <- colnames(inframe.fus)
total <- unique(rbind(inframe.fus[,cols], to.add[,cols], two.callers[,cols], rec.fusions[,cols], rec.fused.genes[,cols]))
total <- total %>% group_by_at((vars(-note))) %>%
  summarise(note = toString(note)) %>%
  as.data.frame() %>%
  unique()

# remove fusions that are in > 1 histology
# not applicable when histology == 1
save(total, file = 'results/fusions_v1.RData')

