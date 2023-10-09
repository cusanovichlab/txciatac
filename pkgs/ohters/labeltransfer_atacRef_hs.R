.libPaths("/home/u23/haozhang1/miniconda2/envs/seurat/lib/R/library")

# from https://github.com/satijalab/seurat/issues/5257
# the error in MapQuery() is caused by an update in uwot
# need to install the previous version of uwot and restart the RStudio session to solve this error.
# devtools::install_version("uwot", version = "0.1.10", repos = "http://cran.us.r-project.org")

#options(stringsAsFactors = FALSE)
library(argparse)
library(tidyverse)
library(Seurat)
library(patchwork)
library(harmony)
library(Signac)
library(S4Vectors)
library(GenomeInfoDb)
#library(EnsDb.Mmusculus.v79)
library(pheatmap)
library(EnsDb.Hsapiens.v86)
#library(rtracklayer)

##################
# # setup variables
# wd_path = "/xdisk/darrenc/haozhang1/scidrop/20210713_scidrop_liver_plus10x/reports"
# setwd(wd_path)
# fragpath = "../fragments/mm10_liver_all.diet_fragments.tsv.gz"
# outfile = "rdata/mm10_liver_allcells_cluster_specific"
# dmat_path = paste0("rdata/mm10_liver_all.diet_cluster.specific_c", 1:6, ".dgcmat.rds")
# peaks = read_tsv("readcount_matrix/mm10_liver_all.diet_cluster.specific_readcounter.matrix.binary.rownames.txt", col_names = F)
# meta_path = readRDS("rdata/mm10_liver_allcells_firstround.seurat.atac.filtered_singlet_harmony.rds")
# # ref_dmat = "/xdisk/darrenc/haozhang1/mouse_atlas_reanalysis/reports/rdata/Liver_62016.dgcmat.rds"
# # ref_meta_path = "/xdisk/darrenc/haozhang1/mouse_atlas_reanalysis/reports/rdata/cell_metadata.txt"
# ref_rna_path = "/xdisk/darrenc/haozhang1/liver_atlas_stst/liver_atlas"
# #ref_rna_meta_path = "/xdisk/darrenc/haozhang1/liver_atlas_stst/annot_mouseStStAll.csv"

parser = argparse::ArgumentParser(description='Script to run first round of clustering .')
# parser$add_argument('--fragpath', required=TRUE, help='Fragment path for reference atac data.')
parser$add_argument('--outfile', required=TRUE, help='Output prefix.')
parser$add_argument('--dmat_path', required=TRUE, help='Prefix of rds file for sparse matrix.')
parser$add_argument('--meta_path', required=TRUE, help='seurat obj from 1st round clustering.')
parser$add_argument('--peak_path', required=TRUE, help='ref peaks of sparse matrix.')
parser$add_argument('--ref_dmat', required=TRUE, help='sparse matrix of reference sparse matrix.')
parser$add_argument('--ref_meta', required=TRUE, help='Reference metedata.')
# parser$add_argument('--tss_path', required=TRUE, help='TSS ref path.')
#parser$add_argument('--hvar', help='Harmony group.by variable.')
#parser$add_argument('--group_var', required=TRUE, default="Rep", help='Group variable for QC violin plot.')
#parser$add_argument('--second_var', help='2nd Variable for UMAP visualization.')
parser$add_argument('--pc', default=30, type="integer", help='Number of PCs used for dimentional reduction.')
parser$add_argument('--res', default=0.8, type="double", help='Resolution for findclusters.')
parser$add_argument('-n', required=TRUE, type="integer", help='Number of rds sparse matrix.')
args = parser$parse_args()

##################
#
# load quary metadata
min_cell = 15
print(paste0("min.cells = ", min_cell))
print(paste0(args$pc, " of PCs used"))
print(paste0(args$res, " of Resolution used"))

##################
# load ref counts
ref_meta <- read_tsv(args$ref_meta)

ref_meta = ref_meta %>%
  column_to_rownames("Barcode")

ref_counts = readRDS(args$ref_dmat)
ref_cells = intersect(rownames(ref_meta), colnames(ref_counts))
ref_counts = ref_counts[, ref_cells]
ref_meta = ref_meta[ref_cells, ]

###########################################
# load quary metadata
metadat = readRDS(args$meta_path)@meta.data
##################
# load quary counts
dmat = paste0(args$dmat_path, 1:args$n, ".dgcmat.rds")
scidrop_count_list = sapply(dmat, readRDS, simplify=FALSE)
scidrop_counts = do.call(cbind, scidrop_count_list)
peaks = read_tsv(args$peak_path, col_names = F)
rownames(scidrop_counts) = peaks$X1

#save sparse matrix
message('Saving read count matrix for query..')
saveRDS(scidrop_counts, file = paste0(args$outfile, "_atacToatac_dgcmat_scidrop.rds"))

# scidrop_counts = readRDS(args$dmat_path)

##############
# Seurat analysis
scipeaks <- StringToGRanges(regions = rownames(scidrop_counts), sep = c("_", "_"))

# check peaks between ref and quary
message('Check peak consistency between ref and quary...')
print(all(rownames(ref_counts) == rownames(scidrop_counts)))

# find peaks that present in at least both 50 ref and quary counts
message('Find high representitive peaks from both ref and quary...')
ref_feature = rownames(ref_counts)[rowSums(ref_counts) >= min_cell]
scidrop_feature = rownames(scidrop_counts)[rowSums(scidrop_counts) >= min_cell]
common_feature = intersect(ref_feature,scidrop_feature)

# subset counts by common_feature
ref_cm_ct = ref_counts[common_feature, ]
scidrop_cm_ct = scidrop_counts[common_feature, ]
print(all(rownames(ref_cm_ct) == rownames(scidrop_cm_ct)))

# make seurat obj
grpeaks <- StringToGRanges(regions = common_feature, sep = c("_", "_"))

message('Total ref cells...')
print(dim(ref_counts))

message('Run Seurat on ref data...')
chrom_ref <- CreateChromatinAssay(
  counts = ref_cm_ct,
  sep = c("_", "-"),
  ranges = grpeaks,
  min.features = 200
)

ref_meta_filter = ref_meta[colnames(chrom_ref), ]

message('Ref cells in seurat obj...')
print(dim(chrom_ref))

seurat_ref <- CreateSeuratObject(
  counts = chrom_ref,
  assay = "peaks",
  meta.data = ref_meta_filter
)

seurat_dim = dim(seurat_ref)
count_dim = dim(ref_counts)
# plot cells filtered
num_cells = rowSums(ref_counts)
num_feature = colSums(ref_counts)
feature_cut = min(rowSums(ref_cm_ct))

message('Feature cut off for Ref cells...')
print(feature_cut)

pdf(paste0(args$outfile, "_atacToatac_cell_sites_histogram_sciATAC_ref.pdf"), 10, 5)
par(mfrow = c(1, 2))
hist(log10(num_cells),main=paste0(c(count_dim[1] - seurat_dim[1]), " Site were removed"), breaks=50)
abline(v=log10(feature_cut),lwd=2,col="indianred")
hist(log10(num_feature),main=paste0(c(count_dim[2] - seurat_dim[2]), " Cells were removed"), breaks=50)
abline(v=log10(200),lwd=2,col="indianred")
dev.off()

# # build annotation
# # extract gene annotations from EnsDb
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# 
# # add the gene information to the object
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg38"
# Annotation(seurat_ref) <- annotations

seurat_ref
seurat_ref[['peaks']]
granges(seurat_ref)

##################################################
message('Clustering reference obj...')

seurat_ref <- FindTopFeatures(seurat_ref, min.cutoff = 'q0')
length(VariableFeatures(seurat_ref))
seurat_ref = RunTFIDF(seurat_ref)
seurat_ref = Signac::RunSVD(seurat_ref)

set.seed(12345)
seurat_ref <- Seurat::RunUMAP(seurat_ref, dims = 2:args$pc, reduction = 'lsi', return.model = TRUE)
seurat_ref <- FindNeighbors(object = seurat_ref, reduction = 'lsi', dims = 2:args$pc)
seurat_ref <- FindClusters(object = seurat_ref, verbose = FALSE, algorithm = 3, resolution = args$res)

message('Run Harmony on reference obj...')
# harmonize lung replicates for bin data
harm.bin <- RunHarmony(
  object = seurat_ref,
  group.by.vars = 'Sample',
  reduction = 'lsi',
  max.iter.cluster = 50,
  max.iter.harmony = 20,
  assay.use = 'peaks',
  project.dim = FALSE
)

set.seed(12345)
harm.bin = Seurat::RunUMAP(harm.bin, dims = 2:args$pc, reduction = 'harmony', return.model = TRUE)
harm.bin = Seurat::FindNeighbors(harm.bin, reduction='harmony', dims = 2:args$pc)
harm.bin = Seurat::FindClusters(harm.bin, verbose = FALSE, algorithm = 3, resolution = args$res)

p_unharm <- DimPlot(seurat_ref, reduction = "umap", group.by = "Sample", pt.size=.1) + ggtitle("bing unharm")

p_harm1 <- DimPlot(harm.bin, reduction = "umap", group.by = "Sample", pt.size=.1) + ggtitle("bing harm")

p_harm2 <- DimPlot(harm.bin, reduction = "umap", pt.size = .1,
                   label = TRUE, repel = TRUE)

p_harm3 <- DimPlot(harm.bin, reduction = "umap", group.by = "Abbreviation", pt.size = .1,
                   label = TRUE, repel = TRUE)

pdf(paste0(args$outfile, "_atacToatac_umap_sciATAC_ref.pdf"), 12, 6)
print(p_unharm | p_harm1)
print(p_harm2 | p_harm3)
dev.off()

message('Saving reference obj...')
saveRDS(harm.bin, file = paste0(args$outfile, '_atacToatac_harmony.seurat.obj_sciATAC_ref.rds'))

#########################################
# process scidrop data
message('Run Seurat on quary data...')
chrom_scidrop <- CreateChromatinAssay(
  counts = scidrop_cm_ct,
  sep = c("_", "-"),
  ranges = grpeaks
)

########################################
message('Check cell consistency between quary matrix and quary metadata...')
print(all(rownames(metadat) == colnames(scidrop_counts)))
print(paste(c("metadata dim = ", dim(metadat)), collapse = " "))
print(paste(c("matrix dim = ", dim(scidrop_counts)), collapse = " "))

if (sum((rownames(metadat) %in% colnames(scidrop_counts))) == dim(metadat)[1]) {
  print("Barcode checking passed.")
} else stop("The barcodes between metadata and count matrix were inconsistent.")

scidrop_meta = metadat[colnames(chrom_scidrop), ]
print(all(rownames(scidrop_meta) == colnames(chrom_scidrop)))

seurat_scidrop <- CreateSeuratObject(
  counts = chrom_scidrop,
  assay = "peaks",
  meta.data = scidrop_meta
)

message('Check quary cell and peak number...')
print(dim(seurat_scidrop))

message('Run SVD on quary data...')
# cell clustering
seurat_scidrop <- FindTopFeatures(seurat_scidrop, min.cutoff = 'q0')
length(VariableFeatures(seurat_scidrop))
seurat_scidrop = RunTFIDF(seurat_scidrop)
seurat_scidrop = Signac::RunSVD(seurat_scidrop)

####################################
#
# label transfer
#
####################################
message('Performing label transfer...')

transfer.anchors <- FindTransferAnchors(
  reference = harm.bin,
  query = seurat_scidrop,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:args$pc
)

message('Saving transfer anchors...')
saveRDS(transfer.anchors, file = paste0(args$outfile, "_atacToatac_transfer.anchors.rds"))

seurat_scidrop <- MapQuery(
  anchorset = transfer.anchors,
  reference = harm.bin,
  query = seurat_scidrop,
  refdata = harm.bin$Abbreviation,
  reference.reduction = "lsi",
  new.reduction.name = "reflsi",
  reduction.model = 'umap'
)

p1 <- DimPlot(harm.bin, reduction = "umap", group.by = "Abbreviation", 
              label = TRUE, repel = TRUE, pt.size=.1) + ggtitle("Reference sciATAC")
p2 <- DimPlot(seurat_scidrop, reduction = "ref.umap", group.by = "predicted.id", 
              label = TRUE, repel = TRUE, pt.size=.1) + NoLegend() + ggtitle("Query sciDrop ATAC")

pdf(paste0(args$outfile, "_atacToatac_projection.umap.pdf"), 12, 6)
print(p1 | p2)
dev.off()

######################################
message('Coembedding reference and quary atac data...')

harm.bin$dataset <- "sciATAC"
seurat_scidrop$dataset <- "sciDrop"
harm.bin$predicted.id = harm.bin$Abbreviation
  
# merge
combined <- merge(harm.bin, seurat_scidrop)
# process the combined dataset
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunTFIDF(combined)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:args$pc)

message("Harmonize ref and quary...")
combined.harm <- RunHarmony(
  object = combined,
  group.by.vars = "dataset",
  reduction = 'lsi',
  max.iter.cluster = 50,
  max.iter.harmony = 20,
  assay.use = 'peaks',
  project.dim = FALSE
)

# re-compute the UMAP using corrected LSI embeddings
set.seed(12345)
combined.harm <- Seurat::RunUMAP(combined.harm, dims = 2:args$pc, reduction = 'harmony')

# create a new UMAP using the integrated embeddings
p1harm <- DimPlot(combined.harm, reduction = "umap", group.by = "dataset") + ggtitle("integrated")
p2harm <- DimPlot(combined, reduction = "umap", group.by = "dataset") + ggtitle("un-integrated")
p3harm <- DimPlot(combined.harm, reduction = "umap", group.by = "predicted.id", split.by = "dataset",
              label = TRUE, repel = TRUE)

pdf(paste0(args$outfile, "_atacToatac_harmony.umap.pdf"), 12, 6)
print(p1harm | p2harm)
print(p3harm)
dev.off()

message('Saving quary obj...')
saveRDS(seurat_scidrop, file = paste0(args$outfile, '_atacToatac_seurat.obj_scidrop.rds'))

message('Saving harmony obj...')
saveRDS(combined.harm, file = paste0(args$outfile, '_atacToatac_seurat.obj_ref_quary_harmony.rds'))

message('Saving query metadata...')
write_tsv(seurat_scidrop@meta.data %>% rownames_to_column("Barcodes"), paste0(args$outfile, '_atacToatac_annotation_metadata_scidrop.txt'))
