.libPaths("/home/u23/haozhang1/miniconda2/envs/seurat/lib/R/library")

# from https://github.com/satijalab/seurat/issues/5257
# the error in MapQuery() is caused by an update in uwot
# need to install the previous version of uwot and restart the RStudio session to solve this error.
# devtools::install_version("uwot", version = "0.1.10", repos = "http://cran.us.r-project.org")

#options(stringsAsFactors = FALSE)
library(argparse)
library(tidyverse)

##################
parser = argparse::ArgumentParser(description='Script to run first round of clustering .')
parser$add_argument('--outfile', required=TRUE, help='Output prefix.')
parser$add_argument('--dmat_path', required=TRUE, help='Prefix of rds file for sparse matrix.')
parser$add_argument('--peak_path', required=TRUE, help='peak file of sparse matrix.')
parser$add_argument('-n', required=TRUE, type="integer", help='Number of rds sparse matrix.')
args = parser$parse_args()

##################
# load quary counts
dmat = paste0(args$dmat_path, 1:args$n, ".dgcmat.rds")
scidrop_count_list = sapply(dmat, readRDS, simplify=FALSE)
scidrop_counts = do.call(cbind, scidrop_count_list)
peaks = read_tsv(args$peak_path, col_names = F)
rownames(scidrop_counts) = peaks$X1

#save sparse matrix
message('Saving read count matrix..')
saveRDS(scidrop_counts, file = paste0(args$outfile, ".dgcmat.rds"))
