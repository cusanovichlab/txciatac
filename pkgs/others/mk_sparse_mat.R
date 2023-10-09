.libPaths("/home/u23/haozhang1/miniconda2/envs/r_env/lib/R/library")
args = commandArgs(TRUE)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(Matrix)

##################
# load binary matrix
bmat = read_tsv(args[1])
#bmat = bmat[, -c(1:3)]

#bmat = bmat %>%
#  column_to_rownames("annot")
#rownames(bmat) = gsub(args[2], "", rownames(bmat))
#print(bmat[1:5, 1:10])
#bmatrix = data.matrix(bmat)
#rm(bmat)
message('Generating sparse matrix...')
counts <- as(data.matrix(bmat), "dgCMatrix")
rm(bmat)
message('Saving sparse matrix...')
saveRDS(counts, file = paste0(args[2], ".dgcmat.rds"))
message('Done saving')
