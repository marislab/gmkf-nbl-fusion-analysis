# Author: Komal S. Rathi
# Date: 08/13/2019
# Function to get pfam domains for fusion genes

setwd('~/Projects/Maris-lab/GMKF_NBL/')

library(biomaRt)
library(dplyr)

# load star fusion results
load('data/star_merged.RData')
star.merge <- as.data.frame(star.merge)
sf <- colsplit(star.merge[,1], '--', names = c("gene1","gene2"))

# load arriba results
load('data/arriba_merged.RData')
ar <- as.data.frame(arriba.merge)
colnames(ar)[1] <- "gene1"
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

# create gene list
genes <- unique(c(ar$gene1, ar$gene2, sf$gene1, sf$gene2))

# get pfam domains for each gene using biomart hg38
grch38 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
la <- listAttributes(grch38)
lf <- listFilters(grch38)
all.hgnc <- getBM(attributes = 'hgnc_symbol', mart = grch38)
genes <- intersect(genes, all.hgnc$hgnc_symbol)
res <- getBM(attributes = c("hgnc_symbol","pfam"),
             filters = "hgnc_symbol", 
             values = genes, 
             mart = grch38)
res.new <- res[res$pfam != "",]
pfam <- read.delim('data/pfamDesc.txt', stringsAsFactors = F)
pfam.res <- merge(res.new, pfam, by.x = 'pfam', by.y = 'pfamAC', all.x = TRUE)

# load('results/pfam_annotated_genes.RData') # all fusion genes annotated using pfam
x <- pfam.res %>% 
  group_by(hgnc_symbol) %>%
  summarise(pfam = toString(unique(pfam)),
            pfamID = toString(unique(pfamID)),
            pfam_desc = toString(unique(pfam_desc))) %>%
  as.data.frame()
x$pfam_domain <- ifelse(grepl('Kinase', x$pfam_desc, ignore.case = TRUE),"pfam_kinase_domain","")
plyr::count(x$pfam_domain)
x$pfam_domain <- ifelse(grepl('transcription', x$pfam_desc, ignore.case = TRUE), "pfam_tf_domain", x$pfam_domain)
plyr::count(x$pfam_domain)
x <- unique(x[,c('hgnc_symbol','pfam_domain')])
pfam.res <- x[x$pfam_domain != "",]
       
save(pfam.res, file = 'results/pfam_annotated_genes.RData')


