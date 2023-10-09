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
#library(S4Vectors)
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
parser$add_argument('--fragpath', required=TRUE, help='Fragment path.')
parser$add_argument('--outfile', required=TRUE, help='Output prefix.')
parser$add_argument('--dmat_path', required=TRUE, help='Prefix of rds file for sparse matrix.')
parser$add_argument('--meta_path', required=TRUE, help='seurat obj from 1st round clustering.')
parser$add_argument('--peak_path', required=TRUE, help='peak file of sparse matrix.')
parser$add_argument('--ref_path', required=TRUE, help='Reference RNA-seq data path.')
parser$add_argument('--hvar', help='Harmony group.by variable.')
parser$add_argument('--group_var', required=TRUE, default="Rep", help='Group variable for QC violin plot.')
parser$add_argument('--second_var', help='2nd Variable for UMAP visualization.')
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

metadat = readRDS(args$meta_path)@meta.data
#rm(meta_path)
metadat = metadat %>%
  dplyr::rename(firstround_cluster = seurat_clusters)

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

# scidrop_counts = readRDS(args$dmat_path)

########################################
#
# since atac atlas liver not good
# use liver atlas rnaseq data to do label transfer
# dataset downloaded from https://www.livercellatlas.org/download.php
#
# use rnaseq_ref.R to pre-process rna seq data
########################################
# Since the cluster not able to handle all cells, only use chow cell to do label transfer
message('Check consistency between metadata and readcount matrix...')
print(all(rownames(metadat) == colnames(scidrop_counts)))
print(paste(c("metadata dim = ", dim(metadat)), collapse = " "))
print(paste(c("matrix dim = ", dim(scidrop_counts)), collapse = " "))

if (sum((rownames(metadat) %in% colnames(scidrop_counts))) == dim(metadat)[1]) {
  print("Barcode checking passed.")
} else stop("The barcodes between metadata and count matrix were inconsistent.")

##############
# Seurat analysis
scipeaks <- StringToGRanges(regions = rownames(scidrop_counts), sep = c("_", "_"))

chrom_assay <- CreateChromatinAssay(
  counts = scidrop_counts,
  sep = c("_", "-"),
  genome = 'hg38',
  ranges = scipeaks,
  fragments = args$fragpath,
  min.cells = min_cell,
  min.features = 200
)

meta_filter = metadat[colnames(chrom_assay), ]

seurat.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = meta_filter
)

message("Cells in Seurat obj...")
print(dim(seurat.atac))
# build annotation
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# add the gene information to the object
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(seurat.atac) <- annotations

seurat.atac
seurat.atac[['peaks']]
granges(seurat.atac)

# cell clustering
seurat.atac = RunTFIDF(seurat.atac)
seurat.atac <- FindTopFeatures(seurat.atac, min.cutoff = 'q0')
print(paste0("Variable length = ", length(VariableFeatures(seurat.atac))))
seurat.atac = Signac::RunSVD(seurat.atac)

set.seed(12345)
seurat.atac = RunUMAP(seurat.atac, reduction = "lsi", dims = 2:args$pc)

########################################################
ref.combined = readRDS(args$ref_path)
ref.combined[['RNA']] <- NULL

#DefaultAssay(ref.combined) <- "integrated"

# calulate gene activity scores
gene.activities <- GeneActivity(seurat.atac, features = VariableFeatures(ref.combined))
# add the gene activity matrix to the Seurat object as a new assay
seurat.atac[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(seurat.atac) <- 'ACTIVITY'

seurat.atac <- NormalizeData(seurat.atac)
seurat.atac <- ScaleData(seurat.atac, features = rownames(seurat.atac))

#saveRDS(seurat.atac, file = paste0(outfile, "_seruat.obj_chow.rds"))

gene_use = VariableFeatures(ref.combined)[which(VariableFeatures(ref.combined) %in% rownames(seurat.atac))]

# Identify anchors
message('Finding anchors..')
transfer.anchors <- FindTransferAnchors(reference = ref.combined, query = seurat.atac, 
                                        features = gene_use,
                                        reference.assay = "integrated", 
                                        query.assay = "ACTIVITY", reduction = "cca")

message('Saving transfer anchors...')
saveRDS(transfer.anchors, file = paste0(args$outfile, "_rnaToatac_transfer.anchors.rds"))

# predict cell type
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = ref.combined$annot,
                                     weight.reduction = seurat.atac[["lsi"]], dims = 2:args$pc)
seurat.atac <- AddMetaData(seurat.atac, metadata = celltype.predictions)

# Co-embedding scRNA-seq and scATAC-seq datasets
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
refdata <- GetAssayData(ref.combined, assay = "integrated", slot = "data")[gene_use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seurat.atac[["lsi"]],
                           dims = 2:args$pc)
seurat.atac[["integrated"]] <- imputation

seurat.atac$annot = seurat.atac$predicted.id
seurat.atac$tech = "scidrop_ATAC"
ref.combined$tech = "scRNA"

coembed <- merge(x = ref.combined, y = seurat.atac)

message("Assay used to combed RNA and ATAC data...")
print(DefaultAssay(coembed))
# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = gene_use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = gene_use, verbose = FALSE)
#coembed <- RunUMAP(coembed, dims = 1:args$pc)

set.seed(2022)
coembed$annot = factor(coembed$annot, levels = sample(unique(coembed$annot)))
coembed$tech = factor(coembed$tech, levels = c("scRNA", "scidrop_ATAC"))

message('Run harmony on combined rna and atac obj...')
coembed.harm <- RunHarmony(
  object = coembed,
  group.by.vars = "tech",
  reduction = 'pca',
  max.iter.cluster = 50,
  max.iter.harmony = 20,
  assay.use = 'integrated',
  project.dim = FALSE
)

set.seed(2022)
coembed.harm <- Seurat::RunUMAP(coembed.harm, dims = 2:args$pc, reduction = 'harmony')

p1pred <- DimPlot(seurat.atac, reduction = "umap", group.by = "predicted.id",
                  label = TRUE, repel = TRUE, pt.size = .1) + NoLegend()
p2pred <- DimPlot(seurat.atac, reduction = "umap", group.by = args$group_var, pt.size=.1)
p3pred <- DimPlot(coembed.harm, group.by = "annot", raster=FALSE,
                  split.by = "tech", label = TRUE, repel = TRUE, pt.size=.1) + NoLegend() + 
  ggtitle("Coembedding RNA and ATAC")

pdf(paste0(args$outfile, "_rnaToatac_umap_coembed.pdf"), 15, 6)
print(p1pred | p2pred)
print(p3pred)
dev.off()

message('Saving RNA and ATAC coembedding obj...')
saveRDS(coembed, file = paste0(args$outfile, "_rnaToatac_coembed_seruat.obj.rds"))

#######################
#
# Clustering scidrop by itself
#
#######################
DefaultAssay(seurat.atac) <- 'peaks'

seurat.atac <- FindNeighbors(object = seurat.atac, reduction = 'lsi', dims = 2:args$pc)
seurat.atac <- FindClusters(object = seurat.atac, verbose = FALSE, algorithm = 3, resolution = args$res)

p1scidrop <- DimPlot(seurat.atac, reduction = "umap", label = TRUE, repel = TRUE, 
                     pt.size = .1)

p2scidrop <- DimPlot(seurat.atac, reduction = "umap", group.by = args$group_var, pt.size=.1)

p3scidrop <- DimPlot(seurat.atac, reduction = "umap", group.by = "predicted.id",
                     label = TRUE, repel = TRUE, pt.size=.1) + NoLegend()

if (!is.null(args$second_var)) {
  p4scidrop <- DimPlot(seurat.atac, reduction = "umap", group.by = args$second_var, pt.size = .1)
  
  pdf(paste0(args$outfile, "_rnaToatac_umap_atac.alone_unintegrated.pdf"), 12, 6)
  print(p1scidrop | p3scidrop)
  print(p2scidrop | p4scidrop)
  dev.off()
} else {
  pdf(paste0(args$outfile, "_rnaToatac_umap_atac.alone_unintegrated.pdf"), 18, 6)
  print(p1scidrop | p2scidrop | p3scidrop)
  dev.off()
}

#################################
#################################
#
# Harmony on Diet and Rep
#
################################
################################
if (!is.null(args$hvar)) {
  message("Running Harmony...")
  
  seurat.atac$unharm_clusters = seurat.atac$seurat_clusters

  print(paste0("Variables used for Harmnoy ", args$hvar))
  print(unique(seurat.atac@meta.data[, args$hvar]))
  
  seurat.harm <- RunHarmony(
    object = seurat.atac,
    group.by.vars = args$hvar,
    reduction = 'lsi',
    max.iter.cluster = 50,
    max.iter.harmony = 20,
    assay.use = 'peaks',
    project.dim = FALSE
  )
  
  depthcorplot = DepthCor(seurat.harm, reduction = 'harmony')
  #elbowplot = ElbowPlot(seurat.harm, reduction = 'harmony', ndims = 50)
  pdf(paste0(args$outfile, ".harmony.depthcor.pdf"), 10, 10)
  print(depthcorplot)
  dev.off()
  
  # re-compute the UMAP using corrected LSI embeddings
  set.seed(12345)
  seurat.harm <- Seurat::RunUMAP(seurat.harm, dims = 2:args$pc, reduction = 'harmony')
  seurat.harm = Seurat::FindNeighbors(seurat.harm, reduction='harmony', dims = 2:args$pc)
  seurat.harm = Seurat::FindClusters(seurat.harm, verbose = FALSE, algorithm = 3, resolution = args$res)
  
  p1harm <- DimPlot(seurat.harm, reduction = "umap", pt.size = .1, 
                    label = TRUE, repel = TRUE) + ggplot2::ggtitle("Harmony integration")
  
  p2harm <- DimPlot(seurat.harm, reduction = "umap", group.by = "unharm_clusters", pt.size = .1,
                    label = TRUE, repel = TRUE) + ggplot2::ggtitle("Unintegrated clusters")
  
  p3harm <- DimPlot(seurat.harm, reduction = "umap", group.by = args$hvar,
                   pt.size=.1)
  
  p4harm <- DimPlot(seurat.harm, reduction = "umap", group.by = "predicted.id", 
                    pt.size = .1, label = TRUE, repel = TRUE)
  
  pdf(paste0(args$outfile, "_rnaToatac_umap_atac.alone_integrated.pdf"), 12, 6)
  print(p1harm | p2harm)
  print(p3harm | p4harm)
  dev.off()
  
  output = paste0(args$outfile, "_rnaToatac_atac.alone_integrated_adj.celltype")
  
  message('Saving unintegrated seurat obj...')
  saveRDS(seurat.atac, file = paste0(args$outfile, "_rnaToatac_atac.alone_unintegrated_seruat.obj.rds"))
  
} else {
  print("No integration performed due to hvar not provided...")
  seurat.harm = seurat.atac
  
  output = paste0(args$outfile, "_rnaToatac_atac.alone_adj.celltype")
}

###########################################################################
###########################################################################

# adjust cell type based on the majority votes for each cluster

###########################################################################
###########################################################################
# use FetchData to pull the identity class (cluster ID), 
# dataset, and original cell labels from novaseq data
# Idents(seurat_harmony) = "seurat_clusters"
message('Adjusting predicted cell types...')
my.data=FetchData(seurat.harm, c("ident", "predicted.id"))

# tabulate cluster ID with cell label
tb = table(my.data$ident, my.data$predicted.id)
tb = as.data.frame.matrix(tb)
#tb = tb %>% dplyr::select(-Unknown)
tb = tb[,order(colSums(tb), decreasing = T)]

pdf(paste0(args$outfile, '_rnaToatac_celltype_confusion.matrix.pdf'), 20, 6)
pheatmap(tb,
         breaks = seq(-2, +2, length = 101),
         angle_col = 45,
         cellwidth = 16,
         cellheight = 8,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "row")
dev.off()

#tb = tb[, !(names(tb) %in% "scidrop")]
# annotate cluster manually with cell label, this aims to annotate cell type for scidrop data
nidents = length(levels(my.data$ident))
celltypes = character(length = nidents)
for (i in 1:nidents) {
  celltypes[i] = colnames(tb)[which.max(tb[i,])]
}

# celltype_idx = apply(tb, 1, which.max)
# #celltype_idx[8] = 4
# #
# pdf(paste0(args$outfile, '.label_transfer.celltype_confusion.matrix.pdf'), 6, 6)
# pheatmap(tb[, unique(celltype_idx)],
#          breaks = seq(-2, +2, length = 101),
#          angle_col = 45,
#          cellwidth = 24,
#          cellheight = 12,
#          cluster_cols = FALSE,
#          cluster_rows = FALSE,
#          scale = "row")
# dev.off()

# add cell type to metadata, since my.data$ident is factor, if ident is 0, as.numeric(my.data$ident) is 1
my.data = rownames_to_column(my.data)
my.data = my.data %>%
  mutate(celltype.adj = celltypes[as.numeric(my.data$ident)])
my.data = my.data %>% column_to_rownames("rowname")
# check if cell type assigned correctly
x = character()
for (i in 1:length(levels(my.data$ident))) {
  x[i] = unique(my.data$celltype.adj[my.data$ident == as.character(i-1)])
}
message('Matching cell types and idents...')
all(x == celltypes)

celltype_meta = my.data %>% dplyr::select(celltype.adj)
# calculate misclassfication rate
#1-mean(my.data2$cell_label == my.data2$celltype_manual)

# add manually assigned cell type to seurat object metadata
if (all(rownames(seurat.harm@meta.data) == rownames(celltype_meta))) {
  seurat.harm = AddMetaData(seurat.harm, celltype_meta)
  print("celltype added")
} else stop("Adjusted cell types were not added because inconsistent rowname.")

if (sum(grepl("\\s", seurat.harm$celltype.adj)) > 0) {
  seurat.harm$celltype.adj = gsub(" ", "_", seurat.harm$celltype.adj)
}

# plot umap with manually assigned cell type and compare to novaseq type
# Idents(seurat_harmony) = "celltype_manual"

p1adj <- DimPlot(seurat.harm, reduction = "umap", 
                  label = TRUE, repel = TRUE, pt.size = .1)

p2adj <- DimPlot(seurat.harm, reduction = "umap", group.by = args$group_var, pt.size = .1)

p3adj <- DimPlot(seurat.harm, reduction = "umap", group.by = "celltype.adj",
                 label = TRUE, repel = TRUE, pt.size = .1) + NoLegend() + 
  ggplot2::ggtitle("Cell type corrected by majority vote")

p4adj <- DimPlot(seurat.harm, reduction = "umap", group.by = "predicted.id", 
              pt.size = .1, label = TRUE, repel = TRUE) + NoLegend() +
  ggplot2::ggtitle("Original predicted cell type")

pdf(paste0(output, "_umap.pdf"), 12, 6)
print(p1adj | p2adj)
print(p3adj | p4adj)
dev.off()

######################################
message("Identify gene markers using gene activity score...")
seurat.harm[['ACTIVITY']] <- NULL

gene.activities <- GeneActivity(seurat.harm)
# add the gene activity matrix to the Seurat object as a new assay
seurat.harm[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(seurat.harm) <- 'ACTIVITY'

seurat.harm <- NormalizeData(
  object = seurat.harm,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat.harm$nCount_ACTIVITY)
)

seurat.harm <- FindVariableFeatures(seurat.harm, selection.method = "vst", nfeatures = 5000)
# scale data with all genes
all.genes <- rownames(seurat.harm)
seurat.harm <- ScaleData(seurat.harm, features = all.genes)

scidrop.markers <- FindAllMarkers(seurat.harm, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
# scidrop.markers %>%
#   group_by(cluster) %>%
#   slice_max(avg_log2FC, n = 10)

scidrop.markers %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  #dplyr::distinct(gene, .keep_all = TRUE) %>%
  group_by(cluster) %>%
  slice_min(p_val_adj, n = 10) %>%
  slice_max(avg_log2FC, n = 10) -> top10

#levels(seurat.harm) <- c(0:24)

top10 = top10 %>%
  arrange(factor(cluster, levels = levels(seurat.harm)))

cluster_heatmap = DoHeatmap(seurat.harm, assay = 'ACTIVITY', 
                            features = top10$gene, 
                            disp.min = -2, disp.max = 2)

pdf(paste0(output, "_cluster.marker.activity_heatmap.pdf"), 16, 10)
print(cluster_heatmap)
dev.off()

# save markers
write_tsv(scidrop.markers %>% rownames_to_column("rowname"), paste0(output, '_marker.genes.txt'))
write_tsv(top10 %>% rownames_to_column("rowname"), paste0(output, '_marker.genes_top10.txt'))

# saveRDS(seurat.harm, file = paste0(args$outfile, '.chow_rna.activity_seurat.obj.rds'))

DefaultAssay(seurat.harm) <- 'peaks'
message("DefaultAssay...")
print(DefaultAssay(seurat.harm))
message("Seurat Idents...")
print(levels(Idents(seurat.harm)))

message('Saving harmony obj...')
saveRDS(seurat.harm, file = paste0(output, '_seurat.obj.rds'))

message('Saving harmony metadata...')
write_tsv(seurat.harm@meta.data %>% rownames_to_column("Barcodes"), paste0(output, '_annotation_metadata.txt'))
