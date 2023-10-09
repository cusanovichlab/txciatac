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
parser$add_argument('--cluster', required=TRUE, help='Cluster of interest.')
parser$add_argument('--group', required=TRUE, help='Column name of metadata containing comparison groups.')
parser$add_argument('--ctl', required=TRUE, help='Control class.')
parser$add_argument('--treat', required=TRUE, help='Treatment class.')
parser$add_argument('--cores', default=1, type="integer", help='Number of cores.')
parser$add_argument('--min_cell', type="integer", help='Number of cores.')
parser$add_argument('--output_directory', required=TRUE, help='Output directory.')
args = parser$parse_args()
#args = list()
#args$monocle_patch = "~/Dropbox/CusanovichLab/Projects/sciDropATAC/monocle_patch.R"
#args$seurat_object = "~/Downloads/liver_mm10_scidrop_seurat.harmony.by.diet_seurat.obj.rds"
#args$cluster = "5"
#args$cores = "1"
#args$output_directory = "~/Dropbox/CusanovichLab/Projects/sciDropATAC/"

source(args$monocle_patch)
seurat_full = readRDS(args$seurat_obj)
message("Subset Seurat obj according to cell bc...")
cells.use = seurat_full@meta.data %>%
  rownames_to_column("bc") %>%
  dplyr::filter(get(args$group) %in% c(args$ctl, args$treat))
message("Check subset groups...")
print("Subset groups...")
print(table(cells.use[, args$group]))
print("All groups...")
print(table(seurat_full@meta.data[,args$group]))

seurat_obj = subset(x = seurat_full, 
                    subset = , cells = cells.use$bc)
seurat_obj_count = GetAssayData(object = seurat_obj, slot = "counts")
clusterfactor = factor((Idents(seurat_obj) == args$cluster)+0,levels=c("0","1"))
message("Observations in cluster of interest...")
print(args$cluster)
print(table(clusterfactor))
print(table(Idents(seurat_obj)))

groupfactor = seurat_obj@meta.data[, args$group]

ctl = seurat_obj_count[,which(clusterfactor=="1" & groupfactor== args$ctl)]
treat = seurat_obj_count[,which(clusterfactor=="1" & groupfactor== args$treat)]
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

#usables = sort(union(which(ctlcount > 5),which(treatcount > 5)))
test_count = cbind(ctl,treat)[usables,]

cells.test = colnames(test_count)

rownames(cells.use) = cells.use$bc
batch_meta = cells.use[colnames(test_count), ]

message("Check batch information...")
print(all(colnames(test_count) == rownames(batch_meta)))

pda = data.frame(CellID = as.character(colnames(test_count)),
                 PeakFrags = as.numeric(log10(colSums(test_count))),
                 Treatment = factor(rep(c(0,1),times=c(ncol(ctl),ncol(treat))),levels=c("0","1")),
                 Batch = as.character(batch_meta$batch))

#names(pda) = c("CellID", "PeakFrags","Treatment")
rownames(pda) = pda$CellID
pda = new("AnnotatedDataFrame", data = pda)

message("Observations in test and control groups...")
print("Treatment...")
print(table(pda$Treatment))
print("Batch...")
print(table(pda$Batch))
message("Confirm observations in test and control groups...")
print(table(Idents(seurat_obj), seurat_obj@meta.data[, args$group]))

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
differtest = differentialGeneTest(submat_cds, fullModelFormulaStr = "~Treatment + Batch + PeakFrags",
                                  reducedModelFormulaStr = "~Batch + PeakFrags", cores=args$cores)

differtest = differtest %>%
  separate(gene_short_name, c("Chr", "Start", "End")) %>%
  arrange(Chr, Start, End) %>%
  rownames_to_column("Peak")

bkg = differtest %>%
  dplyr::select(Chr, Start, End)

message("Differential testing done...")
write_tsv(differtest, paste0(args$output_directory, args$cluster,"_da.by.", args$group, "_all.txt"))
write_tsv(bkg, paste0(args$output_directory, args$cluster,"_da.by.", args$group, "_all_peaks.bed"), col_names = F)

differtest_sig = differtest[which(differtest$qval <= 0.05),]
if(nrow(differtest_sig) > 0){
  differtest_o = differtest_sig[order(differtest_sig$qval),]
  differtest_up = differtest_o[which(differtest_o$beta >0),]
  differtest_down = differtest_o[which(differtest_o$beta <0),]
  #sigers = strsplit2(differtest_o$gene_short_name,"-")
  sigers = differtest_o %>%
    dplyr::select(Chr, Start, End)
  write_tsv(differtest_o, paste0(args$output_directory, args$cluster,"_da.by.", args$group, "_sig.txt"))
  write_tsv(sigers, paste0(args$output_directory, args$cluster,"_da.by.", args$group, "_sig_peaks.bed"), col_names = F)
  if(nrow(differtest_up) > 0){
    #sigers_up = strsplit2(differtest_up$gene_short_name,"-")
    sigers_up = differtest_up %>%
      dplyr::select(Chr, Start, End)
    write_tsv(differtest_up, paste0(args$output_directory, args$cluster,"_da.by.", args$group, "_up.txt"))
    write_tsv(sigers_up, paste0(args$output_directory, args$cluster,"_da.by.", args$group, "_up_peaks.bed"), col_names = F)
  }
  if(nrow(differtest_down) > 0){
    #sigers_down = strsplit2(differtest_down$gene_short_name,"-")
    sigers_down = differtest_down %>%
      dplyr::select(Chr, Start, End)
    write_tsv(differtest_down, paste0(args$output_directory, args$cluster,"_da.by.", args$group, "_down.txt"))
    write_tsv(sigers_down, paste0(args$output_directory, args$cluster,"_da.by.", args$group, "_down_peaks.bed"), col_names = F)
  }
}
