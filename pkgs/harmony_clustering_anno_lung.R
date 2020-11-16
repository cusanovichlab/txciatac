setwd("C:\\Users\\hzhan\\Documents\\Dr.Cusanovich\\data\\sci_atac_drop\\scidrop_100vs800nmSBS")

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

# load mm index
sbs100_index = read_tsv("clustering/scidrop.SBS100_mm_lung_frip0.2.indextable.txt", col_names = F)
sbs800_index = read_tsv("clustering/scidrop.SBS800_mm_lung_frip0.2.indextable.txt", col_names = F)
atlas_index = read_tsv("clustering/cell_metadata.tissue_freq_filtered.txt", col_names = T)

atlas_lung_index = atlas_index %>%
  filter(tissue == "Lung")

###################################################
###################################################
# functions
# reduction means the dimentinal reduction method; small resolution gives you fewer clusters: the seurat tutorial is 0.4 to 1.2
lsi_workflow = function(bmat, dims, metadata=NULL, log_scale_tf=TRUE, reduction='pca.l2', resolution=0.3) {
  tfidf_mat = tfidf(bmat, log_scale_tf=log_scale_tf)
  pca_mat = do_pca(tfidf_mat, dims=max(dims))
  
  seurat_obj = run_dim_reduction(bmat, pca_mat, dims, metadata) %>%
    Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)
  return(seurat_obj)
}

tfidf = function(bmat, frequencies=TRUE, log_scale_tf=TRUE, scale_factor=100000) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  
  # TF-IDF
  tf_idf_counts = safe_tfidf_multiply(tf, idf)
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  return(tf_idf_counts)
}

safe_tfidf_multiply = function(tf, idf) {
  tf = t(tf)
  tf@x <- tf@x * rep.int(idf, diff(tf@p))
  tf = t(tf)
  return(tf)
}

do_pca = function(mat, dims=50) {
  pca.results = irlba(t(mat), nv=dims)
  final_result = pca.results$u %*% diag(pca.results$d)
  rownames(final_result) = colnames(mat)
  colnames(final_result) = paste0('PC_', 1:dims)
  return(final_result)
}

run_dim_reduction = function(atac_matrix, cell_embeddings, dims, metadata=NULL, reduction='pca.l2') {
  if (is.null(metadata)) {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix)
  } else {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix, meta.data = metadata)
  }
  
  seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=cell_embeddings, key='PC_', assay='RNA')
  seurat_obj = seurat_obj %>%
    Seurat::L2Dim(reduction='pca') %>%
    Seurat::RunUMAP(reduction = reduction, dims = dims) %>%
    Seurat::RunTSNE(reduction = reduction, dims = dims) %>%
    Seurat::FindNeighbors(reduction=reduction, nn.eps=0.25, dims=dims)
  return(seurat_obj)
}

# load lung binary matrix
# load atlas data
atlas.bmat_rep1 = read_tsv("clustering/Lung1_62216_mouse.atlas_readcounter.matrix.binary.txt")
atlas.bmat_rep2 = read_tsv("clustering/Lung2_62216_mouse.atlas_readcounter.matrix.binary.txt")

atlas.index = atlas_lung_index %>%
  mutate(treat = tissue.replicate, sample = "atlas_lung") %>%
  select(cell, sample, treat, cell_label)

colnames(sbs100_index) = c("cell", "treat", "sample")
colnames(sbs800_index) = c("cell", "treat", "sample")

sbs100_index = sbs100_index %>%
  mutate(treat = "sbs100", cell_label = "sbs100") %>%
  select(cell, sample, treat, cell_label)

sbs800_index = sbs800_index %>%
  mutate(treat = "sbs800", cell_label = "sbs800") %>%
  select(cell, sample, treat, cell_label)

mm_index_filter = rbind(sbs100_index, sbs800_index, atlas.index)

sbs100_bmat = read_tsv("clustering/scidrop.SBS100.mm_lung_bulk.peaks_readcounter.matrix.binary.txt")

sbs800_bmat = read_tsv("clustering/scidrop.SBS800.mm_lung_bulk.peaks_readcounter.matrix.binary.txt")


atlas.bmat_rep1[1:4, 1:5]
atlas.bmat_rep2[1:4, 1:5]
sbs100_bmat[1:4, 1:5]
sbs800_bmat[1:4, 1:5]

sbs100_bmat$annot = gsub(patter = "mm9", replacement = "",
                          x = sbs100_bmat$annot)
sbs800_bmat$annot = gsub(patter = "mm9", replacement = "",
                         x = sbs800_bmat$annot)

sbs100_bmat[1:2, 1:6]
sbs100.count = sbs100_bmat %>%
  column_to_rownames("annot") %>%
  select(-chr, -start, -end)

sbs800.count = sbs800_bmat %>%
  column_to_rownames("annot") %>%
  select(-chr, -start, -end)

atlas.count_rep1 = atlas.bmat_rep1 %>%
  column_to_rownames("annot") %>%
  select(-chr, -start, -end)

atlas.count_rep2 = atlas.bmat_rep2 %>%
  column_to_rownames("annot") %>%
  select(-chr, -start, -end)
atlas.count <- cbind(atlas.count_rep1, atlas.count_rep2)

sbs100.count_cut = sbs100.count[, mm_index_filter$cell[mm_index_filter$treat == "sbs100"]]
dim(sbs100.count_cut)[2] == sum(mm_index_filter$treat == "sbs100")
sbs100.matrix = data.matrix(sbs100.count_cut)
sbs100.dmat <- as(sbs100.matrix, "dgCMatrix")
rm(sbs100_bmat, sbs100.count, sbs100.count_cut, sbs100.matrix)

sbs800.count_cut = sbs800.count[, mm_index_filter$cell[mm_index_filter$treat == "sbs800"]]
dim(sbs800.count_cut)[2] == sum(mm_index_filter$treat == "sbs800")
sbs800.matrix = data.matrix(sbs800.count_cut)
sbs800.dmat <- as(sbs800.matrix, "dgCMatrix")
rm(sbs800_bmat, sbs800.count, sbs800.count_cut, sbs800.matrix)

atlas.matrix = data.matrix(atlas.count)
atlas.dmat <- as(atlas.matrix, "dgCMatrix")
rm(atlas.bmat_rep1, atlas.bmat_rep2, atlas.count, atlas.matrix, atlas.count_rep1, atlas.count_rep2)

save(atlas.dmat, sbs100.dmat, sbs800.dmat, file = "clustering/lung_sbs_atlas_dgcmat.RData")

# check rowname consistency between atlas and sci drop
sum(row.names(atlas.dmat) != row.names(sbs100.dmat))
sum(row.names(sbs100.dmat) != row.names(sbs800.dmat))

# combine atlas and sci drop
merge.dmat <- cbind(sbs100.dmat, sbs800.dmat, atlas.dmat)
merge.dmat[1:3, 1:3]

# add barcode to rowname of metadata
mm_index_filter = mm_index_filter %>%
  column_to_rownames("cell")

# select top 100000 sites
top_sites = rank(-Matrix::rowSums(merge.dmat)) <= 100000
merge.dmat2 = merge.dmat[top_sites,]
dim(merge.dmat2)

merge_sbs = merge.dmat2[, rownames(mm_index_filter)[which(mm_index_filter$treat %in% c("sbs100", "sbs800"))]]

# save rds files
# save(merge.dmat2, bo.dmat, decoy.dmat, p7.dmat, file = "mm_lung_dgcmatrix.RData")

# calculate term frequency, the default is divide each row by a number
tf = t(t(merge.dmat2) / Matrix::colSums(merge.dmat2))

tf_sbs = t(t(merge_sbs) / Matrix::colSums(merge_sbs))

# compare TF vs log TF
comparison_df = rbind(data.frame(value=tf@x, method='TF'), data.frame(value=log1p(tf@x * 100000), method='log TF'))
ggplot(comparison_df, aes(value)) +
  geom_histogram(bins = 70, aes(fill=method), color='black') +
  theme_classic() +
  scale_fill_manual(values=c("TF"="red", "log TF"="#d3d3d3")) +
  xlab('value in matrix') +
  facet_wrap(~method, scales='free')

comparison_df_sbs = rbind(data.frame(value=tf_sbs@x, method='TF'), data.frame(value=log1p(tf_sbs@x * 100000), method='log TF'))
ggplot(comparison_df_sbs, aes(value)) +
  geom_histogram(bins = 70, aes(fill=method), color='black') +
  theme_classic() +
  scale_fill_manual(values=c("TF"="red", "log TF"="#d3d3d3")) +
  xlab('value in matrix') +
  facet_wrap(~method, scales='free')

# use dim 1 to 50 to make seurat object, due to dim1 is highly affect by sequencing depth
# mouse.seurat.lsi = lsi_workflow(merge.dmat2, dims=1:50, metadata=merge.metdata, log_scale_tf=FALSE, reduction='pca')
mouse.seurat.lsi_log = lsi_workflow(merge.dmat2, dims=2:50, metadata=mm_index_filter, log_scale_tf=TRUE, reduction='pca.l2')

mouse.seurat.lsi_log = SetIdent(mouse.seurat.lsi_log, value="treat")
Idents(mouse.seurat.lsi_log) = factor(Idents(mouse.seurat.lsi_log), levels = c("Lung1_62216", "Lung2_62216", "sbs100", "sbs800"))

p1 = DimPlot(mouse.seurat.lsi_log, reduction = 'umap', pt.size=0.25)

# mouse.seurat.lsi_log_celllabel = SetIdent(mouse.seurat.lsi_log, value="cell_label")
p2 = DimPlot(mouse.seurat.lsi_log, reduction = 'umap', pt.size=0.25, group.by = "cell_label")

pdf("clustering/mouse.lung_logTF_no_integration_umap.pdf", 12, 5)
p1 + p2
dev.off()

mm_index_filter$treat = gsub("100", "", mm_index_filter$treat)
mm_index_filter$treat = gsub("800", "", mm_index_filter$treat)
mm_index_filter$treat = gsub("Lung1_62216", "atlas", mm_index_filter$treat)
mm_index_filter$treat = gsub("Lung2_62216", "atlas", mm_index_filter$treat)

#######################################################################################
#######################################################################################

# run harmony

#######################################################################################
#######################################################################################
tfidf_mat = tfidf(merge.dmat2, log_scale_tf=TRUE)
pca_mat = do_pca(tfidf_mat, dims=50)
pca_mat_l2 = Seurat:::L2Norm(pca_mat[, 2:50])
colnames(pca_mat_l2) = paste0('PCL2_', 1:ncol(pca_mat_l2))

# harmony tutorial use normalized_counts, why we use normalized pc matrix?
sample_count = length(unique(mm_index_filter$treat))
if (sample_count > 1) {
  #message(glue('Running harmony on {sample_count} samples...'))
  pca_embeddings.harmony = harmony::HarmonyMatrix(pca_mat_l2, mm_index_filter$treat, do_pca=FALSE)
} else {
  #message(glue('Only one sample, so not running harmony...'))
  pca_embeddings.harmony = pca_mat_l2
}

colnames(pca_embeddings.harmony) = paste0('PCL2HARMONY_', 1:ncol(pca_embeddings.harmony))
# rownames(pca_embeddings.harmony) = rownames(pca_mat)
all(rownames(pca_embeddings.harmony) == rownames(pca_mat))

seurat_harmony = Seurat::CreateSeuratObject(merge.dmat2, meta.data = mm_index_filter)

seurat_harmony[['pca']] = Seurat::CreateDimReducObject(embeddings=pca_mat, key='PC_', assay='RNA')
seurat_harmony[['pca.l2']] = Seurat::CreateDimReducObject(embeddings=pca_mat_l2, key='PCL2_', assay='RNA')
seurat_harmony[['pca.l2.harmony']] = Seurat::CreateDimReducObject(embeddings=pca_embeddings.harmony, key='PCL2HARMONY_', assay='RNA')

# Running dimensionality reduction and clustering...
seurat_harmony = seurat_harmony %>%
  Seurat::RunUMAP(reduction = 'pca.l2.harmony', dims = 1:49) %>%
  Seurat::RunTSNE(reduction = 'pca.l2.harmony', dims = 1:49) %>%
  Seurat::FindNeighbors(reduction='pca.l2.harmony', nn.eps=0.25, dims=1:49)

## Add UMAP coords to metadata just to make sure
seurat_harmony@meta.data$tissue_umap_1 = Embeddings(seurat_harmony, reduction='umap')[, 1]
seurat_harmony@meta.data$tissue_umap_2 = Embeddings(seurat_harmony, reduction='umap')[, 2]

# parameter used in "run_combined_dim_reduction_LSIwHarmony.R"
# seurat_obj = seurat_obj %>%
#  Seurat::FindClusters(reduction='pca.l2.harmony', n.start=20, resolution=0.3, dims=1:49, min.cluster.size=10)
seurat_harmony = seurat_harmony %>%
  Seurat::FindClusters(reduction='pca.l2.harmony', resolution=0.5, dims=1:49)

###################################################################################
###################################################################################

# plot UMAP

###################################################################################
###################################################################################
# seurat_harmony = SetIdent(seurat_harmony, value="treat")
Idents(seurat_harmony) = "cell_label"
# relevel ident to plot small sample size on top of the plot
# Idents(seurat_harmony) = factor(Idents(seurat_harmony), levels = c("atlas", "sbs"))

p1 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25, split.by = "treat", label = T) + NoLegend()
# Seurat::DimPlot(seurat_harmony, reduction='tsne', pt.size=0.25)

# seurat_harmony_celllabel = SetIdent(seurat_harmony, value="cell_label")
p2 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25, split.by = "sample")

p3 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25, split.by = "treat", 
                     group.by = "seurat_clusters", label = T) + NoLegend()

# harmony
pdf("clustering/mm.lung_logTF.harmony_umap_treat.pdf", 9, 5)
p1
dev.off()

pdf("clustering/mm.lung_logTF.harmony_umap_sample.pdf", 15, 5)
p2
dev.off()

pdf("clustering/mm.lung_logTF.harmony_umap_cluster.pdf", 9, 5)
p3
dev.off()

Seurat::DimPlot(seurat_harmony, reduction='tsne', pt.size=0.25)

save(seurat_harmony, mouse.seurat.lsi_log, file = 'clustering/mm_lung.seurat.RData')

# write.table(seurat_harmony@meta.data, "clustering/mm_liver.metadata.txt", sep='\t', quote = F)

###########################################################################
###########################################################################

# manuall assign cell types

###########################################################################
###########################################################################
# use FetchData to pull the identity class (cluster ID), 
# PC1 scores, dataset, and original cell labels from altas dataset
Idents(seurat_harmony) = "seurat_clusters"
my.data=FetchData(seurat_harmony, c("ident","sample","treat", "cell_label"))
# tabulate cluster ID with cell label
tb = table(my.data$ident, my.data$cell_label)
tb = as.data.frame.matrix(tb)

# annotate cluster manually with cell label, this aims to annotate cell type for scidrop data
celltype = character()
for (i in 1:length(levels(my.data$ident))) {
  celltype[i] = colnames(tb)[tb[i,] == max(tb[i,])]
}

# manual correct celltype
celltype[4] = "Unknown"
celltype[5] = "Unknown"
celltype[7] = "Endothelial I cells"

celltype
# add cell type to metadata, since my.data$ident is factor, if ident is 0, as.numeric(my.data$ident) is 1
my.data = rownames_to_column(my.data)
my.data2 = my.data %>%
  mutate( celltype_manual = celltype[as.numeric(my.data$ident)])

# check if cell type assigned correctly
x = character()
for (i in 1:length(levels(my.data2$ident))) {
  x[i] = unique(my.data2$celltype_manual[my.data2$ident == as.character(i-1)])
}
all(x == celltype)

# calculate misclassfication rate
1-mean(my.data2$cell_label == my.data2$celltype_manual)

# add manually assigned cell type to seurat object metadata
seurat_harmony@meta.data$celltype_manual = my.data2$celltype_manual

# save metadata
write.table(seurat_harmony@meta.data, "clustering/mm_lung.metadata_manual_anno.txt", quote = F, sep = "\t")

# plot umap with manually assigned cell type and compare to altas cell type
Idents(seurat_harmony) = "celltype_manual"
p1 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25,
                     split.by = "treat", label = T) + NoLegend()

p3 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25, split.by = "sample")

p4 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25,
                     group.by = "seurat_clusters", split.by = "treat", label = T) + NoLegend()

pdf("clustering/mm.lung_logTF.harmony_umap_treat_manual_label.pdf", 9, 5)
p1
dev.off()

pdf("clustering/mm.lung_logTF.harmony_umap_sample_manual_label.pdf", 15, 5)
p3
dev.off()

pdf("clustering/mm.lung_logTF.harmony_umap_cluster_manual_label.pdf", 9, 5)
p4
dev.off()

# match cluster color to cell type color
Idents(seurat_harmony) = "seurat_clusters"
names(celltype) = levels(seurat_harmony)

seurat_harmony <- RenameIdents(seurat_harmony, celltype)

Idents(seurat_harmony) = "celltype_manual"
p1 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25, split.by = "treat")
# Seurat::DimPlot(seurat_harmony, reduction='tsne', pt.size=0.25)

# seurat_harmony_celllabel = SetIdent(seurat_harmony, value="cell_label")
p2 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25, split.by = "sample")

Idents(seurat_harmony) = "seurat_clusters"
p3 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25, split.by = "treat")

#####################################################
#####################################################
# check proportion of each cell types
#####################################################
#####################################################

celltype_freq = my.data2 %>%
  filter(treat == "sbs") %>%
  group_by(sample, celltype_manual) %>% 
  tally() %>%
  ungroup() %>%
  select(everything(), count = n) %>%
  group_by(sample) %>%
  dplyr::mutate(tot = sum(count), percentage = count/tot*100)

# reorder cell types based on their percentage in atlas data
# celltype_freq$celltype_manual = factor(celltype_freq$celltype_manual, 
#                                       levels = as.character(celltype_freq$celltype_manual[celltype_freq$treat == "atlas"][order(celltype_freq$percentage[celltype_freq$treat == "atlas"], decreasing = T)]))

celltype_levels = c("WT_lung", "DANZ2_WT_lung","COPD_WT_lung", "HF_lung", "DANZ2_HF_lung", "COPD_HF_lung")
celltype_freq$sample = factor(celltype_freq$sample, levels = celltype_levels)

g = ggplot(celltype_freq, aes(x = celltype_manual, y = percentage, fill=sample))

# use position_dodge to increase the space within groups, 
# and use width to increase the space between groups
pdf("clustering/lung_cell_type_percent.pdf", 8, 6) 
g + geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5) +
  ylab("Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, colour = "black"))
dev.off()
# color by diet

color = rep(c("#F8766D", "#00BFC4"), each = 3)
pdf("clustering/lung_cell_type_percent.pdf", 8, 6) 
g + geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5) +
  ylab("Percentage") +
  scale_fill_manual(values=color) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, colour = "black"))
dev.off()