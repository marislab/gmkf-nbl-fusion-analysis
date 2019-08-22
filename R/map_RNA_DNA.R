# Author: Komal S. Rathi
# Date: 08/20/2019
# Function: Map DNA to RNA sample names

setwd('~/Projects/Maris-lab/GMKF_NBL/')

load('data/clin_mapped.RData')
fus <- read.csv('data/1565033226111-manifest_fusions.csv')
fus <- unique(fus[,c("name","Kids.First.Participant.ID","sample_id","Family.ID")])
colnames(fus) <- paste0('RNA_', colnames(fus))
fus$RNA_name <- gsub('.arriba.fusions.pdf','',fus$RNA_name)


strvar <- read.csv('data/1566319527849-manifest_strvar.csv')
strvar <- unique(strvar[,c("name","Kids.First.Participant.ID","sample_id")])
colnames(strvar) <- paste0('DNA_', colnames(strvar))
strvar$DNA_name <- gsub('.manta_somatic.vep.maf','',strvar$DNA_name)

total <- merge(fus, strvar, by.x = "RNA_Kids.First.Participant.ID", by.y = "DNA_Kids.First.Participant.ID", all.x = TRUE)
to.rm <- setdiff(strvar$DNA_name, total$DNA_name)
for(i in 1:length(to.rm)){
  cmd <- paste0('rm data/manta/',to.rm[i],'*.maf')
  system(cmd)
}
write.table(total, file = 'data/fusions_strvar_mapping.txt', quote = F, sep = "\t", row.names = F)
