setwd("C:\\Users\\hzhan\\Documents\\Dr.Cusanovich\\data\\sci_atac_drop\\test_blocking_conditions\\mouse_lung_clustering")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(Matrix)
library(DelayedArray)
library(Seurat)
library(ggplot2)
library(irlba)
library(patchwork)
library(plyr)
library(stringr)
library(harmony)
library(Signac)

# load dgcmatrix
load("mm_lung_dgcmatrix.RData")

# merge scidrop data
merge.dmat <- cbind(bo.dmat, decoy.dmat, p7.dmat)
merge.dmat[1:3, 1:3]

# load metadata with manual annotation
metadata = read.delim("mm_lung.metadata_manual_anno.txt", sep='\t')

metadata.scidrop = metadata %>%
  filter(treat != "atlas") %>%
  select(cell, treat, celltype_manual)
rm(metadata)

# add cell to rowname
rownames(metadata.scidrop) = metadata.scidrop$cell

# select top 100000 sites
top_sites = rank(-Matrix::rowSums(merge.dmat)) <= 100000
merge.dmat2 = merge.dmat[top_sites,]
dim(merge.dmat2)

# match matrix to metadata
all(rownames(metadata.scidrop) == colnames(merge.dmat2))

# create seurat object
seurat.obj <- CreateSeuratObject(counts = merge.dmat2, meta.data = metadata.scidrop)

# calculate variable features
VariableFeatures(seurat.obj) <- names(which(Matrix::rowSums(seurat.obj) > 11))
length(VariableFeatures(seurat.obj))

# Signac::RunTFIDF, method 3: The log-TF method
seurat.obj = Signac::RunTFIDF(seurat.obj, method = 3)
seurat.obj = Signac::RunSVD(seurat.obj, reduction.name = "svd.m3", reduction.key = "SVDm3_", n = 50)
set.seed(12345)

# find most variable PCs
ElbowPlot(seurat.obj, ndims = 50, reduction = "svd.m3")
# pc = 8, which may not ideal for sciATAC

seurat.obj = Seurat::RunUMAP(seurat.obj, reduction = "svd.m3", dims = 2:20)
seurat.obj = Seurat::RunTSNE(seurat.obj, reduction = "svd.m3", dims = 2:20)
  
# find clusters
seurat.obj <- FindNeighbors(object = seurat.obj, reduction = "svd.m3", dims = 2:20)
seurat.obj <- FindClusters(object = seurat.obj, verbose = FALSE, n.start=20, resolution=.7)

p7 <- DimPlot(seurat.obj, reduction = "umap", group.by = "treat", pt.size=0.5)
p8 <- DimPlot(seurat.obj, reduction = "umap", group.by = "celltype_manual", label = TRUE, repel = TRUE, pt.size=0.5)
p9 <- DimPlot(seurat.obj, reduction = "umap", group.by = "seurat_clusters", pt.size=0.5, label = TRUE)
p10 <- DimPlot(seurat.obj, reduction = "tsne", group.by = "celltype_manual", label = TRUE, repel = TRUE, pt.size=0.5)

pdf("scidrop_mm_lung_umap.pdf", 14, 10)
(p8 + p9) / (p7 + p10)
dev.off()

# split datasets
p11 = DimPlot(seurat.obj, reduction='umap', pt.size=0.5, group.by = "celltype_manual", split.by = "treat")
pdf("scidrop_mm_lung_split_umap.pdf", 12, 5)
p11
dev.off()
