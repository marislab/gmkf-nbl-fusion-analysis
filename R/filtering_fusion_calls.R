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
source('R/merge.cytoband.info.R')
source('R/chrom.plot.R')
source('R/filter.fusAnnot.R')

############## Step 1 (format and merge data) ############## 

# load data
load('data/arriba_merged.RData')
load('data/star_merged.RData')
sf <- as.data.frame(star.merge)
ar <- as.data.frame(arriba.merge)
rm(star.merge, arriba.merge)

# cols to keep
cols <- c('GeneA','GeneB','FusionName','Sample','Caller','Fusion_Type','JunctionReadCount','SpanningFragCount','Confidence','GeneA_bp','GeneB_bp','annots')

# star fusion
# format star fusion output
sf <- format.star(sf)
sf.total <- unique(sf[,cols])

# arriba fusion
# format arriba output for intergenic fusions + add output of fusion annotator to arriba
arriba.FusAnnot <- read.delim('results/arriba_fusions_annotated.txt', stringsAsFactors = F)
ar <- format.arriba(ar, arriba.FusAnnot)
ar.total <- unique(ar[,cols])

# histology (no sub-histology for NBL)
load('data/clin_mapped.RData')

# merge both callers (n = 5062)
all.callers <- rbind(ar.total, sf.total)

# merge callers and clinical information
all.callers <- merge(all.callers, clin, by.x = "Sample", by.y = "name")
all.callers <- unique(all.callers)
colnames(all.callers)[colnames(all.callers) == "FusionName"] <- "Fused_Genes"
rm(ar, sf, ar.total, sf.total, arriba.FusAnnot, clin, cols)

############## Step 1 (format and merge data) ############## 

# Plot of unique fusion calls per chromosome per variable
cytoband <- read.delim('data/cytoBand.txt', stringsAsFactors = F, header = F)
fusioncalls.cyto <- chrom.plot(fusion.calls = all.callers, filepath = 'results', cytoband = cytoband)

############## Step 2 (Apply filters) ############## 

# F1: filter circular RNAs (n = 4846)
all.callers <- all.callers[which(all.callers$GeneA != all.callers$GeneB),]

# F2: filter Red Herring fusions (n = 4716)
# keep fusion relevant to cancer biology where confidence is high
# Using: https://github.com/FusionAnnotator/CTAT_HumanFusionLib/wiki
labels.to.keep <- "Mitelman|chimerdb_omim|chimerdb_pubmed|ChimerKB|ChimerPub|ChimerSeq|Cosmic|YOSHIHARA_TCGA|Klijn_CellLines|Larsson_TCGA|CCLE|HaasMedCancer|GUO2018CR_TCGA|TumorFusionsNAR2018|TCGA_StarF2019|CCLE_StarF2019|TCGA|Oncogene|ArcherDX_panel|FoundationOne_panel|OncocartaV1_panel|OncomapV4_panel"
red.herrings <- "GTEx|GTEx_recurrent|BodyMap|DGD_PARALOGS|HGNC_GENEFAM|Greger_Normal|Babiceanu_Normal|ConjoinG|SELFIE"
f1.red.herrings <- filter.fusAnnot(caller = all.callers, red.herring = red.herrings, labels.to.keep)
all.callers <- all.callers[-which(all.callers$Fused_Genes %in% f1.red.herrings),]

# Filter F3: Remove Read-through fusions (n = 3609)
read.throughs <- "NEIGHBORS|NEIGHBORS_OVERLAP"
read.throughs <- filter.fusAnnot(caller = all.callers, red.herring = read.throughs, labels.to.keep)
all.callers <- all.callers[-which(all.callers$Fused_Genes %in% read.throughs),]

# Filter F0: remove all non-expressed genes
# all.callers <- filter.exp(all.callers)

# Filter F1: remove FusionAnnotator red herrings (n = 2511)
# all.callers <- all.callers[-which(all.callers$Fused_Genes %in% fusions.to.remove),]

# Filter F2: filter Circular RNAs (n = 1850)
# all.callers <- cbind(all.callers, colsplit(all.callers$Fused_Genes, pattern = '--', names = c("GeneA","GeneB")))
# all.callers <- all.callers[!(all.callers$GeneA == all.callers$GeneB & 
#                         all.callers$Confidence != "high" & 
#                         all.callers$Fusion_Type != "in-frame"),]
# circ.rna <- all.callers[which(all.callers$GeneA == all.callers$GeneB),]
# write.table(circ.rna, file = 'results/circ_RNA_txt', quote = F, sep = "\t", row.names = F)
# all.callers <- all.callers[which(all.callers$GeneA != all.callers$GeneB),]

# Filter F3: reduce duplications (n = 3596)
# reduce duplicate fusions to keep one fusion per caller with the highest JunctionReadCount and Confidence
# assign 3 to high confidence, 2 is medium and 1 is low
all.callers$Confidence <- ifelse(all.callers$Confidence %in% "high", 3, ifelse(all.callers$Confidence == "medium", 2, 1))
all.callers <- all.callers %>%
  group_by_at(vars(-JunctionReadCount, -SpanningFragCount, -Confidence)) %>%
  mutate(max.junc = ifelse(JunctionReadCount == max(JunctionReadCount), 'yes', 'no'),
         max.frag = ifelse(SpanningFragCount == max(SpanningFragCount), 'yes', 'no'),
         max.conf = ifelse(Confidence == max(Confidence), 'yes', 'no')) %>%
  filter(max.junc == "yes" & max.frag == "yes" & max.conf == "yes") %>%
  as.data.frame()

# Filter F4: annotate gene1 and gene2 with biotype
# filter any fusion where none of the genes is protein coding
# (n = 3386)
# load('results/fusions_v1.RData')
annot <- read.delim('data/gencode.v27.primary_assembly.annotation.txt', stringsAsFactors = F)
annot <- annot %>% group_by(gene_symbol) %>%
  summarise(gene_id = toString(unique(gene_id)),
            biotype = toString(unique(biotype))) %>%
  unique() %>%
  as.data.frame()
all.callers <- add.biotypes(df = all.callers, annot = annot)

############## Step 2 (Apply filters) ############## 

############## Step 3 (Create lists) ############## 
lists.cols <- c('Fused_Genes','Sample','GeneA_bp','GeneB_bp','Fusion_Type','Caller','note')

# List L1: keep all in-frame fusions
# (n = 174)
inframe.fus <- all.callers %>% 
  filter(Fusion_Type == "in-frame") %>%
  group_by(Fused_Genes, Sample, Fusion_Type, GeneA_bp, GeneB_bp) %>%
  arrange(Caller) %>%
  summarise(Caller = toString(Caller)) %>%
  unique() %>%
  as.data.frame()
inframe.fus$note <- "In-frame fusions"
inframe.fus <- unique(inframe.fus[,lists.cols])
nrow(inframe.fus)
write.table(inframe.fus, file = 'results/inframe_fusions_L1.txt', quote = F, sep = "\t", row.names = F)


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
# (n = 391)
genes.to.search <- unique(c(tsgs, onco))
genes.to.search <- c(paste0('^',genes.to.search,'--'), paste0('--',genes.to.search,'$'))
tsgs_onco <- all.callers[unlist(lapply(genes.to.search, function(x) grep(x, all.callers$Fused_Genes))),]
fusion <- all.callers[unlist(lapply(fusions.to.search, function(x) grep(x, all.callers$Fused_Genes))),]
tsgs_onco_fusion <- unique(rbind(tsgs_onco, fusion))
if(exists('lit.genes.add')){
  to.add <- rbind(lit.genes.add, tsgs_onco_fusion)
} else {
  to.add <- tsgs_onco_fusion
}
to.add <- to.add %>% 
  group_by(Fused_Genes, Sample, Fusion_Type, GeneA_bp, GeneB_bp) %>%
  arrange(Caller) %>%
  summarise(Caller = toString(Caller)) %>%
  unique() %>%
  as.data.frame()
to.add$note <- "Found in onco/tsgs/tcga list"
to.add <- unique(to.add[,lists.cols])
nrow(to.add)
write.table(to.add, file = 'results/Onco_TSG_TCGA_fusions_L2.txt', quote = F, sep = "\t", row.names = F)


# List L3: Found by two callers in the same sample and same breakpoint
# Fusion Type is well annotated
# (n = 169)
two.callers <- all.callers %>%
  group_by(Sample, Fused_Genes, GeneA_bp, GeneB_bp) %>%
  unique() %>%
  arrange(Caller, Fusion_Type) %>%
  summarise(Caller = toString(Caller), 
            Fusion_Type = toString(unique(Fusion_Type)),
            caller.count = n()) %>%
  filter(caller.count >= 2) %>%
  unique() %>%
  as.data.frame()
two.callers$note <- "All callers"
two.callers <- unique(two.callers[,lists.cols])
nrow(two.callers)
write.table(two.callers, file = 'results/BothCallers_fusions_L3.txt', quote = F, sep = "\t", row.names = F)


# List L4: Found in at least 2 samples of the same histology and Fusion Type is well annotated
# For recurrent fusions, filter out low confidence calls
# (n = 227)
rec.fusions <- all.callers %>%
  arrange(Caller, Fusion_Type) %>%
  group_by(Fused_Genes, GeneA_bp, GeneB_bp) %>%
  summarise(Caller = toString(unique(Caller)),
            Fusion_Type = toString(unique(Fusion_Type)),
            Sample = toString(unique(Sample)),
            Sample.count = n()) %>%
  filter(Sample.count > 1) %>%
  unique() %>%
  as.data.frame()
rec.fusions$note <- "Recurrent fusions"
rec.fusions <- unique(rec.fusions[,lists.cols])
nrow(rec.fusions)
rec.fusions <- rec.fusions %>%
  mutate(Sample = strsplit(as.character(Sample), ",")) %>%
  unnest(Sample) %>%
  as.data.frame()
write.table(rec.fusions, file = "results/Recurrent_Fusions_L4.txt", quote = F, sep = "\t", row.names = F)

# List L5: GeneB or GeneA gene recurrently fused within a histology (>= 5 genes)
# For recurrently fused genes, filter out low confidence calls
# (n = 130)
# rec <- cbind(all.callers, colsplit(all.callers$Fused_Genes, pattern = '--', names = c("GeneA","GeneB")))
all.callers$Histology.Broad <- "Neuroblastoma"
rec <- all.callers
# (n = 25)
rec.geneA <- rec %>%
  filter(Confidence > 1) %>%
  group_by(Histology.Broad, Fused_Genes, GeneA, Caller, Fusion_Type) %>%
  dplyr::select(Histology.Broad, Fused_Genes, GeneA, Caller, Fusion_Type) %>%
  unique() %>% 
  group_by(Histology.Broad, GeneA, Caller, Fusion_Type) %>%
  summarise(GeneA.ct = n()) %>%
  filter(GeneA.ct >= 3) %>% as.data.frame()

# (n = 6)
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
rec.fused.genes <- unique(rec.fused.genes[,lists.cols])
nrow(rec.fused.genes)
write.table(rec.fused.genes, file = "results/Recurrently_Fused_Genes_L5.txt", quote = F, sep = "\t", row.names = F)

# merge all five lists (n = 2741)
# collapse by note
total <- unique(rbind(inframe.fus, to.add, two.callers, rec.fusions, rec.fused.genes))
total <- total %>% group_by_at((vars(-note))) %>%
  summarise(note = toString(note)) %>%
  as.data.frame() %>%
  unique()

# remove fusions that are in > 1 histology
# not applicable when histology == 1
save(total, file = 'results/fusions_v1.RData')

