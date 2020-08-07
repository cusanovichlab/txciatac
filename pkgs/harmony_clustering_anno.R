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

# load readcounts
bo_lung_readcounts = read.table("Blocking_oligo_lung.readcounts.report.txt", row.names = 1, sep = '\t')
decoy_lung_readcounts = read.table("DecoyDNA_lung.readcounts.report.txt", row.names = 1, sep = '\t')
p7linear_lung_readcounts = read.table("P7_linear_lung.readcounts.report.txt", row.names = 1, sep = '\t')
p7expo_lung_readcounts = read.table("P7_expo_lung.readcounts.report.txt", row.names = 1, sep = '\t')

p7_lung_readcounts = rbind(p7linear_lung_readcounts, p7expo_lung_readcounts)

# load mm index
bo_mm_index = read_tsv("Blocking_oligo_mm_lung.indextable.txt", col_names = F)
decoy_mm_index = read_tsv("DecoyDNA_mm_lung.indextable.txt", col_names = F)
p7_mm_index = read_tsv("P7_mm_lung.indextable.txt", col_names = F)

# extract mm lung readcounts
bo_mm_readcounts = bo_lung_readcounts[bo_mm_index$X1, ]
decoy_mm_readcounts = decoy_lung_readcounts[decoy_mm_index$X1, ]
p7_mm_readcounts = p7_lung_readcounts[p7_mm_index$X1, ]

rm(bo_lung_readcounts, decoy_lung_readcounts,
   p7linear_lung_readcounts, p7expo_lung_readcounts,
   p7_lung_readcounts)

# caculate frips
bo_mm_readcounts = bo_mm_readcounts %>%
  rownames_to_column() %>%
  mutate(mm_frips = barnyard.mm9.lung.atac_peaks.bed_DHS/mm9) %>%
  mutate(treat = "blocking_oligo")
decoy_mm_readcounts = decoy_mm_readcounts %>%
  rownames_to_column() %>%
  mutate(mm_frips = barnyard.mm9.lung.atac_peaks.bed_DHS/mm9) %>%
  mutate(treat = "decoy")
p7_mm_readcounts = p7_mm_readcounts %>%
  rownames_to_column() %>%
  mutate(mm_frips = barnyard.mm9.lung.atac_peaks.bed_DHS/mm9) %>%
  mutate(treat = "p7")

mm_frip = data.frame(Frips = c(bo_mm_readcounts$mm_frips,
                     decoy_mm_readcounts$mm_frips,
                     p7_mm_readcounts$mm_frips),
                     treat = rep(c("blocking_oligo", "decoy", "p7"), 
                                 times = c(length(bo_mm_readcounts$mm_frips), 
                                           length(decoy_mm_readcounts$mm_frips), 
                                           length(p7_mm_readcounts$mm_frips))))

pdf("Frips_mm_lung.pdf", 6, 8)
ggplot(mm_frip, aes(x = treat, y = Frips)) +
  geom_violin(aes(color = treat), size = 1.5, trim = F) +
  geom_jitter(aes(color = treat), shape=16, position=position_jitter(0.2), alpha = 0.5) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  theme_classic()
dev.off()

# remove frips lower than 0.2
bo_mm_readcounts_cut = bo_mm_readcounts %>%
  filter(mm_frips >= 0.2)
decoy_mm_readcounts_cut = decoy_mm_readcounts %>%
  filter(mm_frips >= 0.2)
p7_mm_readcounts_cut = p7_mm_readcounts %>%
  filter(mm_frips >= 0.2)

mm_readcounts_cut = rbind(bo_mm_readcounts_cut, 
                          decoy_mm_readcounts_cut,
                          p7_mm_readcounts_cut)
mm_index_filter = mm_readcounts_cut %>%
  select(cell = rowname, treat) %>%
  mutate(cell_label = "scidrop")

# load atlas data
atlas.bmat = readRDS("mouse.atlas.qc_filtered_lung.cells.rds")
atlas.metadata = read.delim("mouse.atlas.lung.metadata.txt", sep='\t')

atlas.index = atlas.metadata %>%
  mutate(treat = "atlas") %>%
  select(cell, treat, cell_label)

mm_index_filter = rbind(mm_index_filter, atlas.index)

rm(bo_mm_readcounts_cut,
   decoy_mm_readcounts_cut,
   p7_mm_readcounts_cut,
   bo_mm_readcounts,
   decoy_mm_readcounts,
   p7_mm_readcounts,
   mm_frip,
   bo_mm_index,
   decoy_mm_index,
   p7_mm_index
)

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
bo.bmat = read.table("Blocking_oligo.Mouse_atlas.peaks_readcounter.matrix.binary.txt", header = TRUE)
bo.bmat$annot = gsub(patter = "mm9", replacement = "",
                          x = bo.bmat$annot)
bo.bmat[1:2, 1:6]
bo.count <- bo.bmat[,-(1:4)]
rownames(bo.count) <- bo.bmat[,4]
bo.count_cut = bo.count[, mm_index_filter$cell[mm_index_filter$treat == "blocking_oligo"]]
dim(bo.count_cut)[2] == sum(mm_index_filter$treat == "blocking_oligo")
bo.matrix = data.matrix(bo.count_cut)
bo.dmat <- as(bo.matrix, "dgCMatrix")
rm(bo.bmat, bo.count, bo.count_cut, bo.matrix)

decoy.bmat = read.table("DecoyDNA.Mouse_atlas.peaks_readcounter.matrix.binary.txt", header = TRUE)
decoy.bmat$annot = gsub(patter = "mm9", replacement = "",
                     x = decoy.bmat$annot)
decoy.bmat[1:2, 1:6]
decoy.count <- decoy.bmat[,-(1:4)]
rownames(decoy.count) <- decoy.bmat[,4]
decoy.count_cut = decoy.count[, mm_index_filter$cell[mm_index_filter$treat == "decoy"]]
dim(decoy.count_cut)[2] == sum(mm_index_filter$treat == "decoy")
decoy.matrix = data.matrix(decoy.count_cut)
decoy.dmat <- as(decoy.matrix, "dgCMatrix")
rm(decoy.bmat, decoy.count, decoy.count_cut, decoy.matrix)

p7.bmat = read.table("P7.Mouse_atlas.peaks_readcounter.matrix.binary.txt", header = TRUE)
p7.bmat$annot = gsub(patter = "mm9", replacement = "",
                        x = p7.bmat$annot)
p7.bmat[1:2, 1:6]
p7.count <- p7.bmat[,-(1:4)]
rownames(p7.count) <- p7.bmat[,4]
p7.count_cut = p7.count[, mm_index_filter$cell[mm_index_filter$treat == "p7"]]
dim(p7.count_cut)[2] == sum(mm_index_filter$treat == "p7")
p7.matrix = data.matrix(p7.count_cut)
p7.dmat <- as(p7.matrix, "dgCMatrix")
rm(p7.bmat, p7.count, p7.count_cut, p7.matrix)

# check rowname consistency between atlas and sci drop
sum(row.names(bo.dmat) != row.names(decoy.dmat))
sum(row.names(decoy.dmat) != row.names(p7.dmat))
sum(row.names(p7.dmat) != row.names(atlas.bmat))

# combine atlas and sci drop
merge.dmat <- cbind(bo.dmat, decoy.dmat, p7.dmat, atlas.bmat)
merge.dmat[1:3, 1:3]

# add barcode to rowname of metadata
rownames(mm_index_filter) = mm_index_filter$cell

# select top 100000 sites
top_sites = rank(-Matrix::rowSums(merge.dmat)) <= 100000
merge.dmat2 = merge.dmat[top_sites,]
dim(merge.dmat2)

# save rds files
save(merge.dmat2, bo.dmat, decoy.dmat, p7.dmat, file = "mm_lung_dgcmatrix.RData")

# calculate term frequency, the default is divide each row by a number
tf = t(t(merge.dmat2) / Matrix::colSums(merge.dmat2))

# compare TF vs log TF
comparison_df = rbind(data.frame(value=tf@x, method='TF'), data.frame(value=log1p(tf@x * 100000), method='log TF'))
ggplot(comparison_df, aes(value)) +
  geom_histogram(bins = 70, aes(fill=method), color='black') +
  theme_classic() +
  scale_fill_manual(values=c("TF"="red", "log TF"="#d3d3d3")) +
  xlab('value in matrix') +
  facet_wrap(~method, scales='free')

# use dim 1 to 50 to make seurat object, due to dim1 is highly affect by sequencing depth
# mouse.seurat.lsi = lsi_workflow(merge.dmat2, dims=1:50, metadata=merge.metdata, log_scale_tf=FALSE, reduction='pca')
mouse.seurat.lsi_log = lsi_workflow(merge.dmat2, dims=2:50, metadata=mm_index_filter, log_scale_tf=TRUE, reduction='pca.l2')

mouse.seurat.lsi_log = SetIdent(mouse.seurat.lsi_log, value="treat")
Idents(mouse.seurat.lsi_log) = factor(Idents(mouse.seurat.lsi_log), levels = c("atlas", "p7", "blocking_oligo", "decoy"))

p1 = DimPlot(mouse.seurat.lsi_log, reduction = 'umap', pt.size=0.25)

Idents(mouse.seurat.lsi_log) = "cell_label"
# mouse.seurat.lsi_log_celllabel = SetIdent(mouse.seurat.lsi_log, value="cell_label")
p2 = DimPlot(mouse.seurat.lsi_log, reduction = 'umap', pt.size=0.25)

pdf("mouse.lung_logTF_no_integration_umap.pdf", 12, 5)
p1 + p2
dev.off()

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
  Seurat::FindClusters(reduction='pca.l2.harmony', n.start=20, resolution=0.3, dims=1:49)

save(seurat_harmony, mouse.seurat.lsi_log, file = 'mm_lung.seurat.RData')

###################################################################################
###################################################################################

# plot UMAP

###################################################################################
###################################################################################
seurat_harmony = SetIdent(seurat_harmony, value="treat")
# relevel ident to plot small sample size on top of the plot
Idents(seurat_harmony) = factor(Idents(seurat_harmony), levels = c("atlas", "p7", "blocking_oligo", "decoy"))

p1 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25)
# Seurat::DimPlot(seurat_harmony, reduction='tsne', pt.size=0.25)

Idents(seurat_harmony) = "cell_label"
# seurat_harmony_celllabel = SetIdent(seurat_harmony, value="cell_label")
p2 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25)

# harmony
pdf("mm.lung_logTF.harmony_umap.pdf", 12, 5)
p1 + p2
dev.off()

Seurat::DimPlot(seurat_harmony, reduction='tsne', pt.size=0.25)

write.table(seurat_harmony@meta.data, "mm_lung.metadata.txt", sep='\t', quote = F)

###########################################################################
###########################################################################

# manuall assign cell types

###########################################################################
###########################################################################
# use FetchData to pull the identity class (cluster ID), 
# PC1 scores, dataset, and original cell labels from altas dataset
Idents(seurat_harmony) = "seurat_clusters"
my.data=FetchData(seurat_harmony, c("ident","PCL2HARMONY_1","treat", "cell_label"))
# tabulate cluster ID with cell label
tb = table(my.data$ident, my.data$cell_label)
tb = as.data.frame.matrix(tb)

# annotate cluster manually with cell label, this aims to annotate cell type for scidrop data
celltype = character()
for (i in 1:length(levels(my.data$ident))) {
  celltype[i] = colnames(tb)[tb[i,] == max(tb[i,])]
}

# manual correct celltype
celltype_correct = c("B cells", "Endothelial II cells", "Type I pneumocytes",
                     "Type II pneumocytes", "Type II pneumocytes", "Endothelial I cells",
                     "T cells", "Alveolar macrophages", "Dendritic cells", "Type I pneumocytes",
                     "Hematopoietic progenitors", "NK cells", "Unknown", "Monocytes", "Unknown")

# add cell type to metadata, since my.data$ident is factor, if ident is 0, as.numeric(my.data$ident) is 1
my.data = rownames_to_column(my.data)
my.data2 = my.data %>%
  mutate( celltype_manual = celltype_correct[as.numeric(my.data$ident)])

# check if cell type assigned correctly
x = character()
for (i in 1:length(levels(my.data2$ident))) {
  x[i] = unique(my.data2$celltype_manual[my.data2$ident == as.character(i-1)])
}
all(x == celltype_correct)

# calculate misclassfication rate
1-mean(my.data2$cell_label == my.data2$celltype_manual)

# add manually assigned cell type to seurat object metadata
seurat_harmony@meta.data$celltype_manual = my.data2$celltype_manual

# save metadata
write.table(seurat_harmony@meta.data, "mm_lung.metadata_manual_anno.txt", quote = F, sep = "\t")

# plot umap with manually assigned cell type and compare to altas cell type
Idents(seurat_harmony) = "celltype_manual"
p1 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25)
Idents(seurat_harmony) = "seurat_clusters"
p2 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25)
p1 + p2

# match cluster color to cell type color
Idents(seurat_harmony) = "seurat_clusters"
names(celltype_correct) = levels(seurat_harmony)

seurat_harmony <- RenameIdents(seurat_harmony, celltype_correct)
p3 = Seurat::DimPlot(seurat_harmony, reduction='umap', pt.size=0.25)
pdf("mm.lung_logTF.harmony_manual_celltype_umap.pdf", 12, 5)
p2 + p3
dev.off()

# split datasets
pdf("mm.lung_logTF.harmony_manual_celltype_split_umap.pdf", 16, 5)
DimPlot(seurat_harmony, reduction='umap', pt.size=0.25, split.by = "treat")
dev.off()

#####################################################
#####################################################
# check proportion of each cell types
#####################################################
#####################################################

celltype_freq = my.data2 %>%
  group_by(treat, celltype_manual) %>% 
  tally() %>%
  ungroup() %>%
  select(everything(), count = n) %>%
  group_by(treat) %>%
  dplyr::mutate(tot = sum(count), percentage = count/tot*100)

# reorder cell types based on their percentage in atlas data
celltype_freq$celltype_manual = factor(celltype_freq$celltype_manual, 
                                       levels = as.character(celltype_freq$celltype_manual[celltype_freq$treat == "atlas"][order(celltype_freq$percentage[celltype_freq$treat == "atlas"], decreasing = T)]))

g = ggplot(celltype_freq, aes(x = celltype_manual, y = percentage, fill=treat))

# use position_dodge to increase the space within groups, 
# and use width to increase the space between groups
pdf("cell_type_percent.pdf", 10, 8) 
g + geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5) +
  ylab("Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, colour = "black"))
dev.off()
