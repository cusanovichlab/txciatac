.libPaths("/home/u23/haozhang1/miniconda2/envs/seurat/lib/R/library")
library(monocle)
library(argparse)
library(limma)
library(Seurat)
library(Signac)
library(tidyverse)

parser = argparse::ArgumentParser(description='Script to standard scATAC differential accessibility in monocle.')
parser$add_argument('--monocle_patch', required=TRUE, help='Monocle patch.')
parser$add_argument('--seurat_obj', required=TRUE, help='Seurat Object.')
parser$add_argument('--cluster', type="integer", required=TRUE, help='Cluster of interest.')
parser$add_argument('--cores', type="integer", default=1, help='Number of cores.')
parser$add_argument('--min_cell', type="integer", help='Number of cores.')
parser$add_argument('--cell_count', required=TRUE, help='Data frame show the min cells between sci and scidrop data for each cluster.')
parser$add_argument('--output_directory', required=TRUE, help='Output directory.')
# parser$add_argument('--lung', action="store_true", default=FALSE, help='Mouse lung data?')
args = parser$parse_args()

source(args$monocle_patch)
cell_count = read_tsv(args$cell_count)
seurat_obj = readRDS(args$seurat_obj)
seurat_obj_count = GetAssayData(object = seurat_obj, slot = "counts")
# clusterfactor = factor((Idents(seurat_obj)== args$cluster)+0,levels=c("0","1"))

## subsample 134 cells from each ctl cluster and combine them
# randofactor = seurat_obj@active.ident
# if(args$lung){
#   randofactor[which(randofactor == 16)] = 2
#   randofactor = droplevels(randofactor)
# }
# randos = matrix(,nrow(seurat_obj_count),0)
# set.seed(2021)
# for(i in 0:(length(levels(randofactor))-1)){
# #  if(i == 15){
#   if(i == args$cluster){
#     next
#   }else{
#     subber = seurat_obj_count[,which(seurat_obj@active.ident == i)]
#     subsub = subber[,sample(ncol(subber),134)]
#     randos = cbind(randos,subsub)
#   }
# }

cluster_cells = WhichCells(seurat_obj, idents = args$cluster)
ctl_cells = WhichCells(seurat_obj, idents = args$cluster, invert = TRUE)
cluster_count = cell_count$cluster_min[cell_count$cluster == args$cluster]
ctl_count = cell_count$ctl_min[cell_count$cluster == args$cluster]

message("Subsample cells for clusters...")
sub_cluster_cells = sample(cluster_cells, cluster_count, replace = F)
sub_ctl_cells = sample(ctl_cells, ctl_count, replace = F)

treat = seurat_obj_count[,sub_cluster_cells]
ctl = seurat_obj_count[,sub_ctl_cells]

print(paste0("Number of cells in cluster of interest = ", ncol(treat)))
print(paste0("Number of cells in control cluster = ", ncol(ctl)))

ctlcount = rowSums(ctl)
treatcount = rowSums(treat)
# only test genes that are detected in a 5% of cells in either group
if (!is.na(args$min_cell)) {
  message(paste0("Using the features oberved in ", args$min_cell, " cells for either group..."))
  usables = sort(union(which(ctlcount >= args$min_cell), which(treatcount >= args$min_cell)))
} else {
  message("Using the features oberved in 5% cells for either group...")
  usables = sort(union(which(ctlcount >= 0.05*ncol(ctl)), which(treatcount >= 0.05*ncol(treat))))
}

message("Number of peaks used for test...")
print(length(usables))

test_count = cbind(ctl,treat)[usables,]
pda = data.frame(as.character(colnames(test_count)),as.numeric(log10(colSums(test_count))),factor(rep(c(0,1),times=c(ncol(ctl),ncol(treat))),levels=c("0","1")))
names(pda) = c("CellID", "PeakFrags","Treatment")
rownames(pda) = pda$CellID
pda = new("AnnotatedDataFrame", data = pda)

message("Observations in test and control groups...")
print(table(pda$Treatment))

message("Confirm observations in test and control groups...")
# print(table(Idents(seurat_obj), seurat_obj@meta.data[, args$group]))
print(cell_count[cell_count$cluster == args$cluster,])

fda = data.frame(gene_short_name = as.character(rownames(test_count)))
#fda$gene_short_name = fda$Peak
#names(fda) = "Peak"
rownames(fda) = fda$gene_short_name 
fda = new("AnnotatedDataFrame", data = fda)
seurat_obj_monocle = test_count

#rownames(seurat_obj_monocle) = NULL
#colnames(seurat_obj_monocle) = NULL
submat_cds =  newCellDataSet(seurat_obj_monocle,
                             featureData = fda,
                             phenoData = pda,
                             expressionFamily=binomialff(),
                             lowerDetectionLimit=1)

submat_cds <- detectGenes(submat_cds, min_expr = 0.1)
pData(submat_cds)$Size_Factor = 1

message("Differential testing by group...")
differtest = differentialGeneTest(submat_cds, fullModelFormulaStr = "~Treatment + PeakFrags",
                                  reducedModelFormulaStr = "~PeakFrags", cores=args$cores)

differtest = differtest %>%
  separate(gene_short_name, c("Chr", "Start", "End")) %>%
  arrange(Chr, Start, End) %>%
  rownames_to_column("Peak")

bkg = differtest %>%
  dplyr::select(Chr, Start, End)

message("Differential testing done...")
write_tsv(differtest, paste0(args$output_directory, "_cluster", args$cluster,"_da_all.txt"))
write_tsv(bkg, paste0(args$output_directory, "_cluster", args$cluster,"_da_all_peaks.bed"), col_names = F)

differtest_sig = differtest[which(differtest$qval <= 0.05),]
if(nrow(differtest_sig) > 0){
  differtest_o = differtest_sig[order(differtest_sig$qval),]
  differtest_up = differtest_o[which(differtest_o$beta >0),]
  differtest_down = differtest_o[which(differtest_o$beta <0),]
  #sigers = strsplit2(differtest_o$gene_short_name,"-")
  sigers = differtest_o %>%
    dplyr::select(Chr, Start, End)
  write_tsv(differtest_o, paste0(args$output_directory, "_cluster", args$cluster,"_da_sig.txt"))
  write_tsv(sigers, paste0(args$output_directory, "_cluster", args$cluster,"_da_sig_peaks.bed"), col_names = F)
  if(nrow(differtest_up) > 0){
    #sigers_up = strsplit2(differtest_up$gene_short_name,"-")
    sigers_up = differtest_up %>%
      dplyr::select(Chr, Start, End)
    write_tsv(differtest_up, paste0(args$output_directory, "_cluster", args$cluster,"_da_up.txt"))
    write_tsv(sigers_up, paste0(args$output_directory, "_cluster", args$cluster,"_da_up_peaks.bed"), col_names = F)
  }
  if(nrow(differtest_down) > 0){
    #sigers_down = strsplit2(differtest_down$gene_short_name,"-")
    sigers_down = differtest_down %>%
      dplyr::select(Chr, Start, End)
    write_tsv(differtest_down, paste0(args$output_directory, "_cluster", args$cluster,"_da_down.txt"))
    write_tsv(sigers_down, paste0(args$output_directory, "_cluster", args$cluster,"_da_down_peaks.bed"), col_names = F)
  }
}
