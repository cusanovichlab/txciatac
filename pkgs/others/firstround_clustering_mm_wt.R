.libPaths("/home/u23/haozhang1/miniconda2/envs/seurat/lib/R/library")

options(stringsAsFactors = FALSE)
library(argparse)
library(tidyverse)
library(Matrix)
library(DelayedArray)
library(Seurat)
library(irlba)
library(patchwork)
library(plyr)
library(stringr)
library(harmony)
library(Signac)
library(S4Vectors)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(rtracklayer)
library(mclust)
source("/groups/darrenc/sbin/scidropatac/atac_scrublet_fun.R")
##################
## setup variables
# workdir = "/xdisk/darrenc/haozhang1/scidrop/20210713_scidrop_liver_plus10x/reports"
# setwd(workdir)
# fragpath = "../fragments/mm10_liver_all.diet_fragments.tsv.gz"
# outfile = "rdata/mm10_liver_allcells_firstround"
# dmat_path = "rdata/mm10_liver_all.dgcmat.rds"
# meta_path = "rdata/mm10_liver_allcells.diet_seurat.metadata.txt"
# tss_path = "/groups/darrenc/references/annotations/mm10/mm10.v23.tss.bed.gz"
# pc = 25
# res = 0.6
# hvar = "sample"
# estimated_doublet_rate = 0.031

parser = argparse::ArgumentParser(description='Script to run first round of clustering .')
parser$add_argument('--fragpath', required=TRUE, help='Fragment path.')
parser$add_argument('--outfile', required=TRUE, help='Output prefix.')
parser$add_argument('--dmat_path', required=TRUE, help='rds file for sparse matrix.')
parser$add_argument('--meta_path', required=TRUE, help='metadata path.')
parser$add_argument('--tss_path', required=TRUE, help='TSS ref path.')
parser$add_argument('--hvar', help='Harmony group.by variable.')
parser$add_argument('--vln_group', default="Rep", help='Group variable for QC violin plot.')
parser$add_argument('--lane_var', help='Variable includes lane information, which used to run scrublet separately.')
parser$add_argument('--cluster_var', help='Variable used to isolate cluster specific cells.')
parser$add_argument('--pc', default=30, type="integer", help='Number of PCs used for dimentional reduction.')
parser$add_argument('--res', default=0.8, type="double", help='Resolution for findclusters.')
parser$add_argument('--estimated_doublet_rate', default=0.05, type="double", help='estimated doublet rate.')
args = parser$parse_args()

# mix_levels = c("COPD_WT_lung", "COPD_HF_lung", "COPD_WT_liver", "COPD_HF_liver", "COPD",
#            "DNAZ2_WT_lung", "DNAZ2_HF_lung", "DNAZ2_WT_liver", "DNAZ2_HF_liver", "DNAZ2")
# condition_levels = c("DNAZ2", "COPD")
print(paste0(args$pc, " of PCs used"))
print(paste0(args$res, " of Resolution used"))
# load read counts
counts = readRDS(args$dmat_path)

# load metadata
meta <- read_tsv(args$meta_path)
#colnames(meta) = c("Barcodes", "Diet", "Rep", "Mix")
meta = meta %>%
  column_to_rownames("Barcode")
#meta$Rep = gsub("M7", "M3", meta$Rep)

message("Check bc consistency between meta and raw count matrix")
print(all(rownames(meta) == colnames(counts)))
print(paste(c("metadata dim = ", dim(meta)), collapse = " "))
print(paste(c("matrix dim = ", dim(counts)), collapse = " "))

if (sum((rownames(meta) %in% colnames(counts))) != dim(meta)[1]) {
  stop("The barcodes between metadata and count matrix were inconsistent.")
} else print("Barcode checking passed.")

counts = counts[, rownames(meta)]
print(all(rownames(meta) == colnames(counts)))
message("Raw count matrix dimentions")
count_dim = dim(counts)
print(count_dim)
# [1] 145474  59379

# make seurat obj
peaks_mm10 <- StringToGRanges(regions = rownames(counts), sep = c("_", "_"))

# a bug for ucsc that paralyzes GenomeInfoDb
# Solution: Instead of calling Seqinfo(genome="mm10) which is the default of seurat
# we call SeqinfoForUCSCGenome("mm10")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c("_", "-"),
  genome = 'mm10',
  ranges = peaks_mm10,
  fragments = args$fragpath,
  min.cells = 50,
  min.features = 200
)

meta_filter = meta[colnames(chrom_assay), ]
message("Check bc consistency between filtered meta and Seurat obj")
print(all(rownames(meta_filter) == colnames(chrom_assay)))

seurat.obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = meta_filter
)

seurat_dim = dim(seurat.obj)

# plot cells filtered
num_cells = rowSums(counts)
num_feature = colSums(counts)
#mincell_cut = min(rowSums(seurat.obj))
#minfeature_cut = min(colSums(seurat.obj))

pdf(paste0(args$outfile, ".cell_sites.histogram.pdf"), 10, 5)
par(mfrow = c(1, 2))
hist(log10(num_cells),main=paste0(c(count_dim[1] - seurat_dim[1]), " Site were removed"), breaks=50)
abline(v=log10(50),lwd=2,col="indianred")
hist(log10(num_feature),main=paste0(c(count_dim[2] - seurat_dim[2]), " Cells were removed"), breaks=50)
abline(v=log10(200),lwd=2,col="indianred")
dev.off()

#seurat.obj$mix = factor(seurat.obj$mix, levels = mix_levels)
#seurat.obj$condition = factor(seurat.obj$condition, levels = condition_levels)

# build annotation
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# # change to UCSC style since the data was mapped to hg19
# # seqlevelsStyle stop working due to ucsc bug
# #seqlevelsStyle(annotations) <- 'UCSC'
# # a temporary solution from https://github.com/Bioconductor/GenomeInfoDb/issues/27
# ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
# seqlevels(annotations) <- ucsc.levels
# genome(annotations) <- "mm10"

# add the gene information to the object
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(seurat.obj) <- annotations

seurat.obj
seurat.obj[['peaks']]
granges(seurat.obj)

#############################################################
# Computing QC Metrics
# compute nucleosome signal score per cell
seurat.obj <- NucleosomeSignal(object = seurat.obj)
seurat.obj$nucleosome_group <- ifelse(seurat.obj$nucleosome_signal > 1, 'NS > 1', 'NS < 1')
print(table(seurat.obj$nucleosome_group))
# need to increase the size of the genomic region 
# that's used for creating the fragment histogram
# since there are not many cells in one of the groups
fragfig = FragmentHistogram(object = seurat.obj,
                            region = "chr1-1-150000000")
pdf(paste0(args$outfile, ".fraghist.pdf"), 6, 6)
print(fragfig)
dev.off()

############################################################
# ## build tss reference
# # 
# gene.ranges <- genes(EnsDb.Mmusculus.v79)
# gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
# tss.ranges <- GRanges(
#   seqnames = seqnames(gene.ranges),
#   ranges = IRanges(start = start(gene.ranges), width = 2),
#   strand = strand(gene.ranges)
# )
# seqlevelsStyle(tss.ranges) <- 'UCSC'
# genome(tss.ranges) <- "mm10"
# tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
tss = read_tsv(args$tss_path, col_names = F)
colnames(tss) = c("chr", "start", "end", "gene_name", "score", "strand")
tss_ranges = makeGRangesFromDataFrame(tss, keep.extra.columns=TRUE)
seqlevelsStyle(tss_ranges) <- 'UCSC'
genome(tss_ranges) <- "mm10"

# # calculate TSS score
# need to make sure the chromosome name of fragment file is start with chr
seurat.obj <- TSSEnrichment(object = seurat.obj, tss.positions = tss_ranges)
# seurat.obj$tss_level <- ifelse(seurat.obj$TSS.enrichment > 20, 'High', ifelse(seurat.obj$TSS.enrichment > 2, 'Medium', 'Low'))
# seurat.obj$high.tss <- ifelse(seurat.obj$TSS.enrichment > 2, 'High', 'Low')
# 
# pdf(paste0(args$outfile, ".TSSenrichment.pdf"), 10, 6)
# par(mfrow=c(2,1))
# print(TSSPlot(seurat.obj, group.by = 'high.tss') + NoLegend())
# print(TSSPlot(seurat.obj, group.by = 'tss_level') + NoLegend())
# dev.off()

message("Cells in seurat obj")
print(seurat_dim)
## [1] 131261  56714
message("Cells after removing low Frip")
print(sum(seurat.obj$Frip >= 0.2))
# [1] 55520
message("Cells after removing low TSS")
print(sum(seurat.obj$TSS.enrichment[seurat.obj$Frip >= 0.2] >= 2))
# [1] 55391
message("Cells after removing over Rip")
print(sum(seurat.obj$Rip[seurat.obj$TSS.enrichment >=2 & seurat.obj$Frip >= 0.2] <= 20000))
# [1] 55387

# fraction of reads in peaks
#seurat.obj$pct_reads_in_peaks <- seurat.obj$peak_region_fragments / seurat.obj$passed_filters * 100

# QC plot
# by barnyard
# qcplot1 = VlnPlot(
#   object = seurat.obj,
#   features = c('Frip', 'TSS.enrichment', 'Totalfrags', 'Rip',
#                'Total', 'nucleosome_signal'),
#   pt.size = 0,
#   ncol = 3,
#   group.by = 'Diet',
#   split.by = 'Rep'
# )

col1 <- c("#f3766e", "#11b4e9", "#9488c1")
col2 = c("#faa818", "#41a30d", "#fbdf72", "#367d7d",  
         "#d33502", "#6ebcbc", "#37526d", "#916848",
         "#f5b390", "#342739", "#bed678","#a6d9ee", 
         "#0d74b6", "#60824f", "#725ca5", "#e0598b")
col = c(col1, col2)

p.Frip <- VlnPlot(seurat.obj, features = "Frip",
                  pt.size = 0,
                  group.by = args$vln_group,
                  #split.by = args$vln_group,
                  cols = col) &
  geom_hline(yintercept = 0.2, linetype="dashed")

p.tss <- VlnPlot(seurat.obj, features = "TSS.enrichment",
                 pt.size = 0,
                 group.by = args$vln_group,
                 #split.by = args$vln_group,
                 cols = col) &
  geom_hline(yintercept = 2, linetype="dashed")

p.rip <- VlnPlot(seurat.obj, features = "Rip",
                 pt.size = 0,
                 group.by = args$vln_group,
                 #split.by = args$vln_group,
                 log = T,
                 cols = col) &
  geom_hline(yintercept = 20000, linetype="dashed")

p.comp <- VlnPlot(seurat.obj, features = "Totalfrags",
                  pt.size = 0,
                  group.by = args$vln_group,
                  #split.by = args$vln_group,
                  log = T,
                  cols = col)

p.dedup <- VlnPlot(seurat.obj, features = "Total",
                   pt.size = 0,
                   group.by = args$vln_group,
                   #split.by = args$vln_group,
                   log = T,
                   cols = col)

p.subnuc <- VlnPlot(seurat.obj, features = "nucleosome_signal",
                    pt.size = 0,
                    group.by = args$vln_group,
                    #split.by = args$vln_group,
                    cols = col)
  
pdf(paste0(args$outfile, ".qc_plot.pdf"), 14, 8)
print(wrap_plots(p.Frip, p.tss, p.rip, p.comp, p.dedup, p.subnuc, ncol = 3))
dev.off()

message("Saving seurat obj...")
saveRDS(seurat.obj, file = paste0(args$outfile, ".seurat.obj.rds"))

##############################################
# filter cell by qc
seurat.filter <- subset(
  x = seurat.obj,
  subset = Rip <= 20000 &
    Frip >= 0.2 &
    TSS.enrichment >= 2)

message("Seurat obj dimentions after filtering low quality cells")
print(dim(seurat.filter))
# [1] 131261  55387
# total 3992 cells removed

message("Saving filtered seurat obj...")
saveRDS(seurat.filter, file = paste0(args$outfile, ".seurat.obj.filtered.rds"))

##########################
#
# scrublet: removing cell doublets
#
##########################
#
cell_metadata.df = seurat.filter@meta.data

if (!is.null(args$lane_var)) {
  message('Performing Scrublet by lane splitting...')
  
  lane_var = args$lane_var
  seurat.list <- SplitObject(seurat.filter, split.by = lane_var)
  
  cells_pass_doublet_list = vector(mode = "list", length = length(seurat.list))
  
  for (i in 1:length(seurat.list)) {
    window_matrix = GetAssayData(object = seurat.list[[i]], slot = "counts")
    
    lane = unique(seurat.list[[i]]@meta.data[, lane_var])
    
    print(paste0("processing lane ", lane, " ..."))
    
    # Now run doublet removal using scrublet-like approach
    message('Running doublet detection...')
    scrublet_result = atac_scrublet(window_matrix, estimated_doublet_rate=args$estimated_doublet_rate, fraction_sim_doublets=0.5, dims=2:args$pc)
    
    doublet_scores = scrublet_result$result
    # doublet_score_threshold = quantile(subset(doublet_scores, !simulated_doublet)$doublet_likelihood, 1 - args$estimated_doublet_rate)
    threshold_exp = scrublet_result$threshold_exp
    threshold_mclust = scrublet_result$threshold_mclust
    
    doublet_hist = ggplot(doublet_scores, aes(doublet_likelihood)) +
      geom_histogram(bins=200, aes(y=..density.., fill=simulated_doublet)) +
      geom_vline(xintercept = threshold_exp, linetype='dashed', color='grey') +
      geom_vline(xintercept = threshold_mclust, linetype='dashed', color='red') +
      #geom_text(aes(x=threshold_exp, label="Empirical\n", y=20), colour="blue", angle=90) +
      #geom_text(aes(x=threshold_mclust, label="Mclust\n", y=20), colour="red", angle=90) +
      facet_wrap(~simulated_doublet, scales = "free_y", nrow=2) +
      theme_classic() +
      ggtitle("Doublet score distribution")
    
    doublet_histlog = ggplot(doublet_scores, aes(doublet_likelihood)) +
      geom_histogram(bins=200, aes(y=..density.., fill=simulated_doublet)) +
      geom_vline(xintercept = threshold_exp, linetype='dashed', color='grey') +
      geom_vline(xintercept = threshold_mclust, linetype='dashed', color='red') +
      #geom_text(aes(x=threshold_exp, label="Empirical\n", y=20), colour="blue", angle=90) +
      #geom_text(aes(x=threshold_mclust, label="Mclust\n", y=20), colour="black", angle=90) +
      #scale_y_log10() +
      scale_x_log10() +
      facet_wrap(~simulated_doublet, scales = "free_y", nrow=2) +
      theme_classic() +
      ggtitle("Log (Doublet score) distribution")
    
    message('threshold_exp and threshold_mclust')
    print(c(threshold_exp, threshold_mclust))
    
    pdf(paste0(args$outfile, '_', lane, '.doublet.distribution.pdf'), width = 6, height = 5)
    #par(mfrow=c(2,1))
    print(doublet_hist)
    print(doublet_histlog)
    dev.off()
    
    # UMAP to visulized simulated doublets in simulated data
    scrublet_obj = scrublet_result$seurat_obj
    # try to figure the UMAP space for ambiguous doublets
    p1doublet <- DimPlot(scrublet_obj, reduction = "umap", group.by = "simulated_doublet", pt.size = .01)
    p2doublet <- DimPlot(scrublet_obj, reduction = "umap", group.by = "prediction", pt.size = .01)
    p3doublet <- DimPlot(scrublet_obj, reduction = "umap", group.by = "simulated_doublet", split.by = "prediction", pt.size = .01) +
      ggplot2::ggtitle("split.by = predicted_doublet & color.by = simulated_doublet")
    p4doublet <- DimPlot(scrublet_obj, reduction = "umap", group.by = "prediction", split.by = "simulated_doublet", pt.size = .01) +
      ggplot2::ggtitle("split.by = simulated_doublet & color.by = predicted_doublet")
    
    pdf(paste0(args$outfile, '_', lane, ".doublet.umap.simulated_cells.pdf"), 12, 6)
    print(p1doublet | p2doublet)
    print(p3doublet)
    print(p4doublet)
    dev.off()
    
    message('Saving scrublet results...')
    saveRDS(scrublet_result, file = paste0(args$outfile, '_', lane, ".scrublet_result.rds"))
    
    # calculate how many doublets removed by each cutoff
    tb = table(doublet_scores$simulated_doublet, doublet_scores$prediction)
    tb = as.data.frame.matrix(tb)
    tb = tb %>% rownames_to_column("simulated_doublet")
    
    write_tsv(tb, file = paste0(args$outfile, '_', lane, ".scrublet_confusion_mat.txt"))
    write_tsv(doublet_scores, file = paste0(args$outfile, '_', lane, '.doublet_scores.txt'))
    
    doublet_score_threshold = threshold_mclust
    
    message(paste0("Use doublet likelihood threshold ", round(doublet_score_threshold, 2)))
    # cells_pass_doublet = subset(doublet_scores, !simulated_doublet & !predicted_doublet)$cell
    cells_pass_doublet_list[[i]] = subset(doublet_scores, !simulated_doublet & doublet_likelihood < doublet_score_threshold)
  }
  
  cells_pass_doublet = bind_rows(cells_pass_doublet_list)$cell
} else {
  message('Performing Scrublet without lane splitting...')
  
  window_matrix = GetAssayData(object = seurat.filter, slot = "counts")
  # num_cells_ncounted = rowSums(window_matrix)
  # ncounts = window_matrix[num_cells_ncounted >= 50,]
  # new_counts = colSums(ncounts)
  # quantile(new_counts,probs=0.01)
  # ncounts = ncounts[,new_counts >= quantile(new_counts,probs=0.01)]
  
  #print(all(rownames(cell_metadata.df) == colnames(window_matrix)))
  
  # Now run doublet removal using scrublet-like approach
  message('Running doublet detection...')
  scrublet_result = atac_scrublet(window_matrix, estimated_doublet_rate=args$estimated_doublet_rate, fraction_sim_doublets=0.5, dims=2:args$pc)
  
  doublet_scores = scrublet_result$result
  # doublet_score_threshold = quantile(subset(doublet_scores, !simulated_doublet)$doublet_likelihood, 1 - args$estimated_doublet_rate)
  threshold_exp = scrublet_result$threshold_exp
  threshold_mclust = scrublet_result$threshold_mclust
  
  doublet_hist = ggplot(doublet_scores, aes(doublet_likelihood)) +
    geom_histogram(bins=200, aes(y=..density.., fill=simulated_doublet)) +
    geom_vline(xintercept = threshold_exp, linetype='dashed', color='grey') +
    geom_vline(xintercept = threshold_mclust, linetype='dashed', color='red') +
    #geom_text(aes(x=threshold_exp, label="Empirical\n", y=20), colour="blue", angle=90) +
    #geom_text(aes(x=threshold_mclust, label="Mclust\n", y=20), colour="red", angle=90) +
    facet_wrap(~simulated_doublet, scales = "free_y", nrow=2) +
    theme_classic() +
    ggtitle("Doublet score distribution")
  
  doublet_histlog = ggplot(doublet_scores, aes(doublet_likelihood)) +
    geom_histogram(bins=200, aes(y=..density.., fill=simulated_doublet)) +
    geom_vline(xintercept = threshold_exp, linetype='dashed', color='grey') +
    geom_vline(xintercept = threshold_mclust, linetype='dashed', color='red') +
    #geom_text(aes(x=threshold_exp, label="Empirical\n", y=20), colour="blue", angle=90) +
    #geom_text(aes(x=threshold_mclust, label="Mclust\n", y=20), colour="black", angle=90) +
    #scale_y_log10() +
    scale_x_log10() +
    facet_wrap(~simulated_doublet, scales = "free_y", nrow=2) +
    theme_classic() +
    ggtitle("Log (Doublet score) distribution")
  
  message('threshold_exp and threshold_mclust')
  print(c(threshold_exp, threshold_mclust))
  ## dim 30
  # 96.9%
  # 0.2083888 0.3026378
  ## dim 20
  # 96.9%
  # 0.1726704 0.2143743
  
  pdf(paste0(args$outfile, '.doublet.distribution.pdf'), width = 6, height = 5)
  #par(mfrow=c(2,1))
  print(doublet_hist)
  print(doublet_histlog)
  dev.off()
  
  # UMAP to visulized simulated doublets in simulated data
  scrublet_obj = scrublet_result$seurat_obj
  # try to figure the UMAP space for ambiguous doublets
  p1doublet <- DimPlot(scrublet_obj, reduction = "umap", group.by = "simulated_doublet", pt.size = .01)
  p2doublet <- DimPlot(scrublet_obj, reduction = "umap", group.by = "prediction", pt.size = .01)
  p3doublet <- DimPlot(scrublet_obj, reduction = "umap", group.by = "simulated_doublet", split.by = "prediction", pt.size = .01) +
    ggplot2::ggtitle("split.by = predicted_doublet & color.by = simulated_doublet")
  p4doublet <- DimPlot(scrublet_obj, reduction = "umap", group.by = "prediction", split.by = "simulated_doublet", pt.size = .01) +
    ggplot2::ggtitle("split.by = simulated_doublet & color.by = predicted_doublet")
  
  pdf(paste0(args$outfile, ".doublet.umap.simulated_cells.pdf"), 12, 6)
  print(p1doublet | p2doublet)
  print(p3doublet)
  print(p4doublet)
  dev.off()
  
  message('Saving scrublet results...')
  saveRDS(scrublet_result, file = paste0(args$outfile, ".scrublet_result.rds"))
  
  # calculate how many doublets removed by each k
  
  tb = table(doublet_scores$simulated_doublet, doublet_scores$prediction)
  tb = as.data.frame.matrix(tb)
  tb = tb %>% rownames_to_column("simulated_doublet")
  
  write_tsv(tb, file = paste0(args$outfile, ".scrublet_confusion_mat.txt"))
  
  # determine the threshold based on the UMAP of simulated doublets
  doublet_score_threshold = threshold_mclust
  
  message(paste0("Use doublet likelihood threshold ", round(doublet_score_threshold, 2)))
  # cells_pass_doublet = subset(doublet_scores, !simulated_doublet & !predicted_doublet)$cell
  cells_pass_doublet = subset(doublet_scores, !simulated_doublet & doublet_likelihood < doublet_score_threshold)$cell
  
  write_tsv(doublet_scores, file = paste0(args$outfile, '.doublet_scores.txt'))
}

doublet_meta = cell_metadata.df %>%
  rownames_to_column("cell") %>%
  mutate(doublet = ifelse(cell %in% cells_pass_doublet, "singlet", "doublet")) %>%
  column_to_rownames("cell")

write_tsv(doublet_meta %>% rownames_to_column("cell"), file = paste0(args$outfile, '.doublet_metadata.txt'))

# UMAP to visualize predicted doublet in read data
doublet_meta = doublet_meta[colnames(seurat.filter), ]
# seurat.filter$doublet = factor(doublet_meta$doublet, levels = c("singlet", "doublet"))
add_doublet = doublet_meta %>% dplyr::select(doublet)
print(all(rownames(add_doublet) == colnames(seurat.filter)))
seurat.filter = AddMetaData(seurat.filter, add_doublet)
  
message("Doublets identified...")
print(table(seurat.filter$doublet))

seurat.filter = RunTFIDF(seurat.filter)
seurat.filter <- FindTopFeatures(seurat.filter, min.cutoff = 'q0')
print(length(VariableFeatures(seurat.filter)))

seurat.filter = Signac::RunSVD(seurat.filter)

set.seed(12345)
seurat.filter = RunUMAP(seurat.filter, reduction = "lsi", dims = 2:args$pc)
#seurat.filter <- RunTSNE(seurat.filter, reduction = "lsi", dims = 2:30)
seurat.filter <- FindNeighbors(object = seurat.filter, reduction = 'lsi', dims = 2:args$pc)
seurat.filter <- FindClusters(object = seurat.filter, algorithm = 3, resolution = args$res)

p1doublet <- DimPlot(seurat.filter, reduction = "umap", group.by = "doublet", pt.size = .1)

p2doublet <- DimPlot(seurat.filter, reduction = "umap", group.by = "seurat_clusters", pt.size=.1)

p3doublet <- DimPlot(seurat.filter, reduction = "umap", group.by = "seurat_clusters",
                     split.by = "doublet", label = TRUE, repel = TRUE, pt.size = .1) + NoLegend()

# p4doublet <- DimPlot(seurat.filter, reduction = "umap", group.by = "doublet", 
#                      split.by = "Diet", label = TRUE, repel = TRUE, pt.size=.1) + NoLegend()

p4doublet <- DimPlot(seurat.filter, reduction = "umap", group.by = "doublet",
                     split.by = args$vln_group, label = TRUE, repel = TRUE, pt.size = .1) + NoLegend()

if (!is.null(args$lane_var)) {
  p5doublet <- DimPlot(seurat.filter, reduction = "umap", group.by = "doublet",
                       split.by = args$lane_var, label = TRUE, repel = TRUE, pt.size = .1) + NoLegend()
 
  pdf(paste0(args$outfile, ".doublet.umap.pdf"), 12, 6)
  print(p1doublet | p2doublet)
  print(p3doublet)
  print(p4doublet)
  print(p5doublet)
  dev.off()
} else {
  pdf(paste0(args$outfile, ".doublet.umap.pdf"), 12, 6)
  print(p1doublet | p2doublet)
  print(p3doublet)
  print(p4doublet)
  dev.off()
}

#############################
# remove doublets
message("Removing doublets...")
seurat.singlet <- subset(
  x = seurat.filter,
  subset = doublet == "singlet")
print(table(seurat.singlet$doublet))
## dim30
# singlet 
#  54415
## dim20
# singlet
# 54152
message("Checking singlet purity...")
print(all(rownames(doublet_meta[doublet_meta$doublet == "singlet", ]) == colnames(seurat.singlet)))
#seurat.singlet$depth <- ifelse(seurat.singlet$peak_region_fragments > 1000, 'High_depth', 'Low_depth')

# re-cluster singlet cells
seurat.singlet = RunTFIDF(seurat.singlet)
seurat.singlet <- FindTopFeatures(seurat.singlet, min.cutoff = 50)
length(VariableFeatures(seurat.singlet))
seurat.singlet = Signac::RunSVD(seurat.singlet)

set.seed(12345)
seurat.singlet = RunUMAP(seurat.singlet, reduction = "lsi", dims = 2:args$pc)
seurat.singlet <- FindNeighbors(object = seurat.singlet, reduction = 'lsi', dims = 2:args$pc)
seurat.singlet <- FindClusters(object = seurat.singlet, algorithm = 3, resolution = args$res)

message("Saving Seurat obj after filtering doublets...")
saveRDS(seurat.singlet, file = paste0(args$outfile, '.seurat.obj.filtered_singlet.rds'))
rm(window_matrix)

p1umap <- DimPlot(seurat.singlet, reduction = "umap", pt.size = .1,
                  label = TRUE, repel = TRUE)

p2umap <- DimPlot(seurat.singlet, reduction = "umap", group.by = args$vln_group, pt.size = .1)

p3umap <- DimPlot(seurat.singlet, reduction = "umap",
                  split.by = args$vln_group, label = TRUE, repel = TRUE, pt.size=.1)

if ("Mix" %in% colnames(meta)) {
  p4umap <- DimPlot(seurat.singlet, reduction = "umap",
                    split.by = "Mix", label = TRUE, repel = TRUE, pt.size=.1)
  
  if (!is.null(args$lane_var)) {
    p5umap <- DimPlot(seurat.singlet, reduction = "umap",
                      group.by = args$lane_var, label = TRUE, repel = TRUE, pt.size=.1)
    
    pdf(paste0(args$outfile, ".singlet_umap.pdf"), 12, 6)
    print(p1umap | p2umap)
    print(p1umap | p5umap)
    print(p3umap)
    print(p4umap)
    dev.off()
  } else {
    pdf(paste0(args$outfile, ".singlet_umap.pdf"), 12, 6)
    print(p1umap | p2umap)
    print(p3umap)
    print(p4umap)
    dev.off()
  }
} else {
  
  if (!is.null(args$lane_var)) {
    p4umap <- DimPlot(seurat.singlet, reduction = "umap",
                      group.by = args$lane_var, label = TRUE, repel = TRUE, pt.size=.1)
    
    pdf(paste0(args$outfile, ".singlet_umap.pdf"), 12, 6)
    print(p1umap | p2umap)
    print(p1umap | p4umap)
    print(p3umap)
    dev.off()
  } else {
    pdf(paste0(args$outfile, ".singlet_umap.pdf"), 12, 6)
    print(p1umap | p2umap)
    print(p3umap)
    dev.off()
  }
}

###############################
#
# Harmony on assay
#
##############################
#dims.use = 2:50,
if (!is.null(args$hvar)) {
  message("Running Harmony...")
  print(all(Idents(seurat.singlet) == seurat.singlet$seurat_clusters))
  seurat.singlet$unharm_clusters = seurat.singlet$seurat_clusters
  
  print(paste0("Variables used for Harmnoy ", args$hvar))
  print(unique(seurat.singlet@meta.data[, args$hvar]))
  
  seurat.harm <- RunHarmony(
    object = seurat.singlet,
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
  seurat.harm = Seurat::FindClusters(seurat.harm, algorithm = 3, resolution = args$res)
  
  p1harm <- DimPlot(seurat.harm, reduction = "umap", pt.size = .1,
                    label = TRUE, repel = TRUE) + ggplot2::ggtitle("Harmony integration")
  
  p2harm <- DimPlot(seurat.harm, reduction = "umap", group.by = args$vln_group, pt.size = .1)
  
  p3harm <- DimPlot(seurat.harm, reduction = "umap",
                    split.by = args$vln_group, label = TRUE, repel = TRUE, pt.size=.1)
  
  p4harm <- DimPlot(seurat.harm, reduction = "umap", pt.size = .1,
                    label = TRUE, repel = TRUE) + NoLegend() + ggplot2::ggtitle("Harmony clusters")
  
  p5harm <- DimPlot(seurat.harm, reduction = "umap", group.by = "unharm_clusters", pt.size = .1,
                    label = TRUE, repel = TRUE) + NoLegend() + ggplot2::ggtitle("Unintegrated clusters")
  
  pdf(paste0(args$outfile, ".harmony.umap.pdf"), 12, 6)
  print(p1harm | p2harm)
  print(p3harm)
  print(p4harm | p5harm)
  dev.off()
  
  message("Saving Harmony Seurat obj...")
  saveRDS(seurat.harm, file = paste0(args$outfile, '.seurat.obj.filtered_singlet_harmony.rds'))
  
} else {
  print("No integration performed.")
  seurat.harm = seurat.singlet
}

if (!is.null(args$cluster_var)) {
  message("Isolating cell bc used for calling cluster and treatment specific peaks...")
  
  cluster_var = args$cluster_var
  seurat.harm$cluster_treat = paste0(seurat.harm@meta.data[,cluster_var],
                                     "_c",
                                     seurat.harm$seurat_clusters)
  cluster_level = levels(factor(seurat.harm$cluster_treat))
  cells = vector(mode = "list", length = length(cluster_level))
  names(cells) = cluster_level
  
  for (i in 1:length(cells)) {
    idx = which(seurat.harm$cluster_treat == names(cells)[i])
    cells[[i]] = data.frame(cells = names(idx),
                            clusters = names(cells)[i])
  }
  
  cell_count = as.data.frame(sapply(cells, nrow)) %>%
    rownames_to_column("Cluster")
  colnames(cell_count) = c('Cluster', "Count")
  
  message('check cell number extracted from each cluster')
  print(all(cell_count$Count == table(seurat.harm$cluster_treat)))
  
  write_tsv(cell_count, paste0(args$outfile, ".cluster_cell.count.txt"))
  
  # merge clusters
  cluster_cells = bind_rows(cells)
  
  write_tsv(cluster_cells, paste0(args$outfile, ".cluster_specific_cells.indextable.txt"), col_names = F)
  
} else{
  message("Isolating cell bc used for calling cluster specific peaks...")
  
  cells = vector(mode = "list", length = length(levels(Idents(seurat.harm))))
  names(cells) = levels(Idents(seurat.harm))
  
  for (i in 1:length(cells)) {
    cells[[i]] = data.frame(cells = WhichCells(seurat.harm, idents = i-1),
                            clusters = paste0("cluster", i-1))
  }
  
  cell_count = as.data.frame(sapply(cells, nrow)) %>%
    rownames_to_column("Cluster")
  colnames(cell_count) = c('Cluster', "Count")
  
  message('check cell number extracted from each cluster')
  print(all(cell_count$Count == table(Idents(seurat.harm))))
  
  write_tsv(cell_count, paste0(args$outfile, ".cluster_cell.count.txt"))
  
  # merge clusters
  cluster_cells = bind_rows(cells)
  
  write_tsv(cluster_cells, paste0(args$outfile, ".cluster_specific_cells.indextable.txt"), col_names = F)
}
