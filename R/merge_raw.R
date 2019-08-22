# Author: Komal S Rathi
# Date: 07/19/2019
# Function to merge expression and fusion data 

setwd('~/Projects/Maris-lab/GMKF_NBL/')

# source libraries
library(R.utils)
library(data.table)
library(reshape2)

# universal function to merge expr and fusion results
merge.res <- function(nm){
  print(nm)
  sample_name <- gsub('.*[/]|[.].*','',nm)
  x <- data.table::fread(nm)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    x$sample_name <- sample_name
    return(x)
  } 
}

# expression (TPM)
dir <- 'data/expr/'
lf <- list.files(path = dir, pattern = "*.gz", recursive = T, full.names = T)
expr.mat <- lapply(lf, FUN = function(x) merge.res(x))
expr.mat <- data.table::rbindlist(expr.mat)
expr.mat <- dcast(expr.mat, gene_id~sample_name, value.var = 'TPM')
expr.mat <- as.data.frame(expr.mat)
expr.mat <- cbind(colsplit(expr.mat$gene_id, '_', names = c("gene_id", "gene_symbol")), expr.mat)
expr.mat <- expr.mat[,-3]

# expression (raw counts)
dir <- 'data/expr/'
lf <- list.files(path = dir, pattern = "*.gz", recursive = T, full.names = T)
expr.mat <- lapply(lf, FUN = function(x) merge.res(x))
expr.mat <- data.table::rbindlist(expr.mat)
expr.mat <- dcast(expr.mat, gene_id~sample_name, value.var = 'expected_count')
expr.mat <- as.data.frame(expr.mat)
expr.mat <- cbind(colsplit(expr.mat$gene_id, '_', names = c("gene_id", "gene_symbol")), expr.mat)
expr.mat <- expr.mat[,-3]
save(expr.mat, file = 'data/expr_counts.RData')

# arriba
dir <- 'data/arriba/'
lf <- list.files(path = dir, pattern = "*.tsv", recursive = T, full.names = T)
arriba.merge <- lapply(lf, FUN = function(x) merge.res(x))
arriba.merge <- data.table::rbindlist(arriba.merge)

# star fusion
dir <- 'data/star/'
lf <- list.files(path = dir, pattern = "*.tsv", recursive = T, full.names = T)
star.merge <- lapply(lf, FUN = function(x) merge.res(x))
star.merge <- data.table::rbindlist(star.merge)

# manta
dir <- 'data/manta/'
lf <- list.files(path = dir, pattern = "*.maf", recursive = T, full.names = T)
manta.merge <- lapply(lf, FUN = function(x) merge.res(x))
manta.merge <- data.table::rbindlist(manta.merge)

# save summarized data
save(expr.mat, file = 'data/expr.RData')
save(arriba.merge, file = 'data/arriba_merged.RData')
save(star.merge, file = 'data/star_merged.RData')
save(manta.merge, file = 'data/manta_merged.RData')
