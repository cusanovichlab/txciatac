#######
# setup functions
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

plot_pc_correlation = function(seurat_obj, reduction, column='nCount_RNA') {
  coords = Seurat::Embeddings(seurat_obj, reduction=reduction)
  column_value = seurat_obj@meta.data[, column]
  correlations = apply(coords, 2, function(x) {cor(x, column_value, method='spearman')})
  correlations_df = data.frame(correlation=correlations, PC=1:ncol(coords))
  
  plot_obj = ggplot(correlations_df, aes(PC, correlation)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 0, linetype='dashed', color='red')
  
  return(plot_obj)
}

plot_clustering_comparison = function(seurat_obj1, seurat_obj2, reduction, description1='', description2 = '', cluster_column1='RNA_snn_res.0.3', cluster_column2='RNA_snn_res.0.3') {
  # Clusters as called on each dataset
  seurat_obj1 = SetIdent(seurat_obj1, value=cluster_column1)
  seurat_obj2 = SetIdent(seurat_obj2, value=cluster_column2)
  
  p1 = DimPlot(seurat_obj1, reduction = 'tsne', pt.size=0.25) +
    ggtitle(description1)
  
  p2 = DimPlot(seurat_obj2, reduction = 'tsne', pt.size=0.25) +
    ggtitle(description2)
  
  # Now swap the labels to verify they are finding the same groups
  seurat_obj1@meta.data$cluster_seurat_obj2 = seurat_obj2@meta.data[, cluster_column2]
  seurat_obj2@meta.data$cluster_seurat_obj1 = seurat_obj1@meta.data[, cluster_column1]
  
  seurat_obj1 = SetIdent(seurat_obj1, value='cluster_seurat_obj2')
  seurat_obj2 = SetIdent(seurat_obj2, value='cluster_seurat_obj1')
  
  p3 = DimPlot(seurat_obj1, reduction = reduction, pt.size=0.25) +
    ggtitle(paste0(description1, ' (', description2, ' clusters)'))
  
  p4 = DimPlot(seurat_obj2, reduction = reduction, pt.size=0.25) +
    ggtitle(paste0(description2, ' (', description1, ' clusters)'))
  
  (p1 + p3) / (p2 + p4)
}

### scrublet from Andrew

# atac_scrublet = function(bmat, k=NULL, fraction_sim_doublets=2, estimated_doublet_rate=0.1, dims=2:30, scale.factor=10000) {
#   # Determine KNN parameters
#   if (is.null(k)) {
#     k = round(0.5 * sqrt(ncol(bmat)))
#   }
#   kadj = round(k * (1+fraction_sim_doublets))
#   
#   # Perform TFIDF on original dataset
#   message('[scrublet atac] Performing LSI-logTF on dataset...')
#   
#   # TF-IDF on original dataset
#   bmat_colsums = Matrix::colSums(bmat)
#   tf_original = safe_column_scale(bmat, bmat_colsums)
#   tf_original@x = log1p(tf_original@x * scale.factor)
#   idf_original = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
#   
#   bmat.tfidf = safe_tfidf_multiply(tf_original, idf_original)
#   rm(tf_original)
#   
#   bmat.pca = sparse_prcomp_irlba(t(bmat.tfidf), n=max(dims), center=FALSE, scale.=FALSE)
#   rm(bmat.tfidf)
#   
#   # Make simulated doublets
#   message('[scrublet atac] Simulating doublets...')
#   set.seed(2019)
#   doublet_half1 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
#   doublet_half2 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
#   doublets.bmat = bmat[, doublet_half1] + bmat[, doublet_half2]
#   doublets.bmat@x[doublets.bmat@x > 1] = 1
#   colnames(doublets.bmat) = paste0('DOUBLET_', 1:ncol(doublets.bmat))
#   
#   # Perform TF-IDF on doublets using IDF from original dataset
#   doublet_colsums = bmat_colsums[doublet_half1] + bmat_colsums[doublet_half2] ## approximate (because of binarization after sum) to save time, but could recalculate
#   tf_doublets = safe_column_scale(doublets.bmat, doublet_colsums)
#   rm(doublets.bmat)
#   tf_doublets@x = log1p(tf_doublets@x * scale.factor)
#   doublets.tfidf = safe_tfidf_multiply(tf_doublets, idf_original)
#   rm(tf_doublets)
#   
#   # Project doublets into PCA space and weight by variance explained
#   message('[scrublet atac] Projecting doublets into PCA space...')
#   doublets.pca = t(doublets.tfidf) %*% bmat.pca$rotation
#   rm(doublets.tfidf)
#   doublets.pca.weighted = doublets.pca
#   
#   # Jam it all into a Seurat object for some of the downstream stepst
#   message('[scrublet atac] Making Seurat object...')
#   rownames(bmat.pca$x) = colnames(bmat)
#   combined.metadata = rbind(data.frame(cell=colnames(bmat), doublet=FALSE), data.frame(cell=rownames(doublets.pca), doublet=TRUE))
#   rownames(combined.metadata) = combined.metadata$cell
#   
#   combined.pca = rbind(bmat.pca$x, doublets.pca.weighted)
#   
#   rownames(combined.pca) = rownames(combined.metadata)
#   combined.seurat = CreateSeuratObject(counts=t(combined.pca), meta.data = combined.metadata)
#   combined.seurat[['pca']] = Seurat::CreateDimReducObject(embeddings=as.matrix(combined.pca), key='PC_', assay='RNA')
#   
#   message('[scrublet atac] Finding KNN...')
#   combined.seurat = combined.seurat %>%
#     Seurat::L2Dim(reduction='pca') %>%
#     Seurat::FindNeighbors(reduction='pca.l2', dims=dims, k=kadj, compute.SNN = FALSE) # nn.eps = 0.25    
#   
#   # From KNN, calculate doublet likelihood as defined in Scrublet paper
#   message('[scrublet atac] Calculating doublet neighbors...')
#   doublet_mask = ifelse(combined.seurat@meta.data$doublet, 1, 0)
#   doublet_neighbors = Matrix::rowSums(safe_column_multiply(combined.seurat@meta.data$RNA_nn, doublet_mask))
#   
#   message('[scrublet atac] Finalizing doublet likelihoods...')
#   doublet_score = doublet_neighbors / kadj
#   q = (doublet_neighbors + 1)/ (kadj + 2)
#   doublet_likelihood = q * estimated_doublet_rate / fraction_sim_doublets / (1 - estimated_doublet_rate - q * (1 - estimated_doublet_rate - estimated_doublet_rate / fraction_sim_doublets))
#   
#   # Return Seurat object with doublet likelihood as an extra column
#   result = data.frame(cell=rownames(combined.seurat@meta.data), 'doublet_score'=doublet_score, 'doublet_likelihood'=doublet_likelihood, 'simulated_doublet'=combined.seurat@meta.data$doublet)
#   rownames(result) = result$cell
#   return(result)
# }

### scrublet using signac TF-IDF

# atac_scrublet = function(bmat, k=NULL, fraction_sim_doublets=2, estimated_doublet_rate=0.1, dims=2:30, scale.factor=10000) {
#   # Determine KNN parameters
#   if (is.null(k)) {
#     k = round(0.5 * sqrt(ncol(bmat)))
#   }
#   kadj = round(k * (1+fraction_sim_doublets))
#   
#   # Perform TFIDF on original dataset
#   message('[scrublet atac] Performing LSI-logTF on dataset...')
#   
#   # TF-IDF on original dataset
#   bmat_colsums = Matrix::colSums(bmat)
#   tf_original = safe_column_scale(bmat, bmat_colsums)
#   idf_original = ncol(bmat) / Matrix::rowSums(bmat)
#   bmat.tfidf = safe_tfidf_multiply(tf_original, idf_original)
#   bmat.tfidf@x = log1p(bmat.tfidf@x * scale.factor)
#   rm(tf_original)
#   
#   bmat.pca = sparse_prcomp_irlba(t(bmat.tfidf), n=max(dims), center=FALSE, scale.=FALSE)
#   rm(bmat.tfidf)
#   
#   # Make simulated doublets
#   message('[scrublet atac] Simulating doublets...')
#   set.seed(2021)
#   doublet_half1 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
#   doublet_half2 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
#   doublets.bmat = bmat[, doublet_half1] + bmat[, doublet_half2]
#   doublets.bmat@x[doublets.bmat@x > 1] = 1
#   colnames(doublets.bmat) = paste0('DOUBLET_', 1:ncol(doublets.bmat))
#   
#   # Perform TF-IDF on doublets using IDF from original dataset
#   doublet_colsums = bmat_colsums[doublet_half1] + bmat_colsums[doublet_half2] ## approximate (because of binarization after sum) to save time, but could recalculate
#   tf_doublets = safe_column_scale(doublets.bmat, doublet_colsums)
#   rm(doublets.bmat)
#   doublets.tfidf = safe_tfidf_multiply(tf_doublets, idf_original)
#   doublets.tfidf@x = log1p(doublets.tfidf@x * scale.factor)
#   rm(tf_doublets)
#   
#   # Project doublets into PCA space and weight by variance explained
#   message('[scrublet atac] Projecting doublets into PCA space...')
#   doublets.pca = t(doublets.tfidf) %*% bmat.pca$rotation
#   rm(doublets.tfidf)
#   doublets.pca.weighted = doublets.pca
#   
#   # Jam it all into a Seurat object for some of the downstream stepst
#   message('[scrublet atac] Making Seurat object...')
#   rownames(bmat.pca$x) = colnames(bmat)
#   combined.metadata = rbind(data.frame(cell=colnames(bmat), doublet=FALSE), data.frame(cell=rownames(doublets.pca), doublet=TRUE))
#   rownames(combined.metadata) = combined.metadata$cell
#   
#   combined.pca = rbind(bmat.pca$x, doublets.pca.weighted)
#   
#   rownames(combined.pca) = rownames(combined.metadata)
#   combined.seurat = CreateSeuratObject(counts=as(t(combined.pca), "dgCMatrix"), meta.data = combined.metadata)
#   combined.seurat[['pca']] = Seurat::CreateDimReducObject(embeddings=as.matrix(combined.pca), key='PC_', assay='RNA')
#   
#   message('[scrublet atac] Finding KNN...')
#   combined.seurat = combined.seurat %>%
#     Seurat::L2Dim(reduction='pca') %>%
#     Seurat::FindNeighbors(reduction='pca.l2', dims=dims, k.param=kadj, compute.SNN = FALSE) # nn.eps = 0.25
#   
#   # From KNN, calculate doublet likelihood as defined in Scrublet paper
#   message('[scrublet atac] Calculating doublet neighbors...')
#   doublet_mask = ifelse(combined.seurat@meta.data$doublet, 1, 0)
#   doublet_neighbors = Matrix::rowSums(safe_column_multiply(combined.seurat@meta.data$RNA_nn, doublet_mask))
#   
#   message('[scrublet atac] Finalizing doublet likelihoods...')
#   doublet_score = doublet_neighbors / kadj
#   q = (doublet_neighbors + 1)/ (kadj + 2)
#   doublet_likelihood = q * estimated_doublet_rate / fraction_sim_doublets / (1 - estimated_doublet_rate - q * (1 - estimated_doublet_rate - estimated_doublet_rate / fraction_sim_doublets))
#   
#   # Return Seurat object with doublet likelihood as an extra column
#   result = data.frame(cell=rownames(combined.seurat@meta.data), 'doublet_score'=doublet_score, 'doublet_likelihood'=doublet_likelihood, 'simulated_doublet'=combined.seurat@meta.data$doublet)
#   rownames(result) = result$cell
#   
#   threshold_exp = quantile(subset(result, !simulated_doublet)$doublet_likelihood, 1 - estimated_doublet_rate)
#   
#   message('[scrublet atac] Estimating doublet threshold using Mclust...')
#   calldouble = Mclust(data.frame(log10(result$doublet_likelihood)),G=2)
#   threshold_mclust = min(result[which(calldouble$classification == 2 & calldouble$uncertainty < 0.05),"doublet_likelihood"])
#   
#   result = result %>%
#     mutate(predicted_exp = ifelse(doublet_likelihood <= threshold_exp, FALSE, TRUE),
#            predicted_mclust = ifelse(doublet_likelihood <= threshold_mclust, FALSE, TRUE)) %>%
#     mutate(prediction = ifelse(!predicted_exp & !predicted_mclust, "singlet", 
#                                ifelse(predicted_exp & predicted_mclust, "doublets", 
#                                       ifelse(!predicted_exp & predicted_mclust, "mclust_doublets", "exp_doublets"))))
#   # do umap on combined pca spaces
#   message('[scrublet atac] Run UMAP...')
#   #print(all(rownames(result) == colnames(combined.seurat)))
#   combined.seurat = AddMetaData(combined.seurat, result)
#   combined.seurat = RunUMAP(combined.seurat, reduction='pca.l2', dims=dims)
#   
#   return(list(result = result, seurat_obj = combined.seurat, threshold_exp = threshold_exp, threshold_mclust = threshold_mclust))
# }

### scrublet using signac TF-IDF and andrew dimentional reduction

atac_scrublet = function(bmat, k=NULL, fraction_sim_doublets=2, estimated_doublet_rate=0.1, dims=2:30, scale.factor=10000) {
  # Determine KNN parameters
  if (is.null(k)) {
    k = round(0.5 * sqrt(ncol(bmat)))
  }
  kadj = round(k * (1+fraction_sim_doublets))
  
  # Perform TFIDF on original dataset
  message('[scrublet atac] Performing LSI-logTF on dataset...')
  
  # TF-IDF on original dataset
  bmat_colsums = Matrix::colSums(bmat)
  tf_original = safe_column_scale(bmat, bmat_colsums)
  idf_original = ncol(bmat) / Matrix::rowSums(bmat)
  bmat.tfidf = safe_tfidf_multiply(tf_original, idf_original)
  bmat.tfidf@x = log1p(bmat.tfidf@x * scale.factor)
  rm(tf_original)
  
  bmat.pca = sparse_prcomp_irlba(t(bmat.tfidf), n=max(dims), center=FALSE, scale.=FALSE)
  rm(bmat.tfidf)
  
  # Make simulated doublets
  message('[scrublet atac] Simulating doublets...')
  set.seed(2021)
  doublet_half1 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
  doublet_half2 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
  doublets.bmat = bmat[, doublet_half1] + bmat[, doublet_half2]
  doublets.bmat@x[doublets.bmat@x > 1] = 1
  colnames(doublets.bmat) = paste0('DOUBLET_', 1:ncol(doublets.bmat))
  
  # Perform TF-IDF on doublets using IDF from original dataset
  doublet_colsums = bmat_colsums[doublet_half1] + bmat_colsums[doublet_half2] ## approximate (because of binarization after sum) to save time, but could recalculate
  tf_doublets = safe_column_scale(doublets.bmat, doublet_colsums)
  rm(doublets.bmat)
  doublets.tfidf = safe_tfidf_multiply(tf_doublets, idf_original)
  doublets.tfidf@x = log1p(doublets.tfidf@x * scale.factor)
  rm(tf_doublets)
  
  # Project doublets into PCA space and weight by variance explained
  message('[scrublet atac] Projecting doublets into PCA space...')
  doublets.pca = t(doublets.tfidf) %*% bmat.pca$rotation
  rm(doublets.tfidf)
  doublets.pca.weighted = doublets.pca
  
  # Jam it all into a Seurat object for some of the downstream stepst
  message('[scrublet atac] Making Seurat object...')
  rownames(bmat.pca$x) = colnames(bmat)
  combined.metadata = rbind(data.frame(cell=colnames(bmat), doublet=FALSE), data.frame(cell=rownames(doublets.pca), doublet=TRUE))
  rownames(combined.metadata) = combined.metadata$cell
  
  combined.pca = rbind(bmat.pca$x, doublets.pca.weighted)
  
  rownames(combined.pca) = rownames(combined.metadata)
  combined.seurat = CreateSeuratObject(counts=as(t(combined.pca), "dgCMatrix"), meta.data = combined.metadata)
  combined.seurat[['pca']] = Seurat::CreateDimReducObject(embeddings=as.matrix(combined.pca), key='PC_', assay='RNA')
  
  message('[scrublet atac] Finding KNN...')
  combined.seurat = combined.seurat %>%
    Seurat::L2Dim(reduction='pca') %>%
    Seurat::FindNeighbors(reduction='pca.l2', dims=dims, k.param=kadj, compute.SNN = FALSE) # nn.eps = 0.25
  
  # From KNN, calculate doublet likelihood as defined in Scrublet paper
  message('[scrublet atac] Calculating doublet neighbors...')
  doublet_mask = ifelse(combined.seurat@meta.data$doublet, 1, 0)
  doublet_neighbors = Matrix::rowSums(safe_column_multiply(combined.seurat@meta.data$RNA_nn, doublet_mask))
  
  message('[scrublet atac] Finalizing doublet likelihoods...')
  doublet_score = doublet_neighbors / kadj
  q = (doublet_neighbors + 1)/ (kadj + 2)
  doublet_likelihood = q * estimated_doublet_rate / fraction_sim_doublets / (1 - estimated_doublet_rate - q * (1 - estimated_doublet_rate - estimated_doublet_rate / fraction_sim_doublets))
  
  # Return Seurat object with doublet likelihood as an extra column
  result = data.frame(cell=rownames(combined.seurat@meta.data), 'doublet_score'=doublet_score, 'doublet_likelihood'=doublet_likelihood, 'simulated_doublet'=combined.seurat@meta.data$doublet)
  rownames(result) = result$cell
  
  # isolated simulated doublets
  sim_doublet = result %>%
    dplyr::filter(simulated_doublet)
  
  threshold_exp = quantile(subset(result, !simulated_doublet)$doublet_likelihood, 1 - estimated_doublet_rate)
  
  message('[scrublet atac] Estimating doublet threshold using Mclust...')
  calldouble = Mclust(data.frame(sim_doublet$doublet_likelihood),G=2)
  threshold_mclust = min(sim_doublet[which(calldouble$classification == 2 & calldouble$uncertainty < 0.05),"doublet_likelihood"])
  
  result = result %>%
    mutate(predicted_exp = ifelse(doublet_likelihood < threshold_exp, FALSE, TRUE),
           predicted_mclust = ifelse(doublet_likelihood < threshold_mclust, FALSE, TRUE)) %>%
    mutate(prediction = ifelse(!predicted_exp & !predicted_mclust, "singlet",
                               ifelse(predicted_exp & predicted_mclust, "common_doublets",
                                      ifelse(!predicted_exp & predicted_mclust, "mclust_doublets", "exp_doublets"))))
  # do umap on combined pca spaces
  message('[scrublet atac] Run UMAP...')
  #print(all(rownames(result) == colnames(combined.seurat)))
  combined.seurat = AddMetaData(combined.seurat, result)
  set.seed(2021)
  combined.seurat = RunUMAP(combined.seurat, reduction='pca.l2', dims=dims)
  
  return(list(result = result, seurat_obj = combined.seurat, threshold_exp = threshold_exp, threshold_mclust = threshold_mclust))
}

safe_column_scale = function(bmat, scale_vector) {
  bmat@x <- bmat@x / rep.int(scale_vector, diff(bmat@p))
  return(bmat)
}

safe_column_multiply = function(bmat, scale_vector) {
  bmat@x <- bmat@x * rep.int(scale_vector, diff(bmat@p))
  return(bmat)
}

sparse_prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, ...) {
  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  
  args <- list(A=x, nv=n, scale.=scale., center=center)
  if (!missing(...)) args <- c(args, list(...))
  s <- do.call(irlba, args=args)
  
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  ans$center <- args$center
  if (retx)
  {
    ans <- c(ans, list(d=s$d, x = s$u %*% diag(s$d)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}
