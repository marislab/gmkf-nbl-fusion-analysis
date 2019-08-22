# Author: Komal S. Rathi
# Date: 08/20/2019
# Function: Map clinical file to RNA sample names

setwd('~/Projects/Maris-lab/GMKF_NBL/')

# clinical file mapped to sample names
clin <- read.csv('data/target_update_gmkf_only_export_6-22-2017.csv')
clin <- unique(clin[,c('usi','inss_stage','mycn_status','risk','diagnosis')])
map <- read.csv('data/1565033226111-manifest.csv')
map <- unique(map[,c("name","Family.ID")])
map$name <- gsub('[.].*', '', map$name)
clin <- merge(clin, map, by.x = 'usi', by.y = 'Family.ID')
save(clin, file = 'data/clin_mapped.RData')
