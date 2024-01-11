#' Collection of functions to assist in analyses
#' @author carlga ~ Carlos Gallardo
#'


#' Performs principal component analysis
#' @param x A dataframe or matrix
#' @param top_var Number of top variable features to use
#' @param center Center data
#' @param scale. Scale data
#' @return List with PCA results
#' 
#' @author carlga
#' 
doPCA <- function(x, top_var = 500, center = T, scale. = F) {
  require(matrixStats)
  
  rvars <- matrixStats::rowVars(as.matrix(x))
  select <- order(rvars, decreasing = TRUE)[seq_len(min(top_var, length(rvars)))]
  pca <- prcomp(t(x[select, ]), center = center, scale. = scale.)
  
  res <- list()
  res[['perc_var']] <- pca$sdev^2/sum(pca$sdev^2)*100
  res[['pcs']] <- as.data.frame(pca$x)
  res[['rotation']] <- as.data.frame(pca$rotation)
  
  return(res)
}


#' Scales data row-wise using selected method
#' @param x A dataframe or matrix
#' @param method Method used for scaling
#' @return Data frame with scaled data
#' 
#' @author carlga
#'
scaleData <- function(x, method = NULL) {
  
  if(is.null(method)) stop('Please, indicate a valid method\n')
  
  if(method == "zscore") {
    row.mean <- apply(x, 1, mean)
    row.sd <- apply(x, 1, sd)
    res <- (x - row.mean) / row.sd
  }
  else if(method == "minmax") {
    row.min <- apply(x, 1, min)
    row.max <- apply(x, 1, max)
    res <- (x - row.min) / (row.max - row.min)
  }
  else if(method == "max") {
    row.max <- apply(x, 1, max)
    res <- x / row.max
  }
  else stop('This method is not available\n')
  
  return(as.data.frame(res))
}


#' Normalizes data column-wise using selected method
#' @param x A dataframe or matrix
#' @param len Vector with gene/transcript lengths (for RPKM/TPM normalization)
#' @param method Method used for normalization
#' @return Data frame with normalized data
#' 
#' @author carlga
#'
normalizeData <- function(x, len = NULL, method = NULL) {
  
  if(is.null(method)) stop('Please, indicate a valid method\n')
  
  x <- as.matrix(x)
  
  if(method == "CPM") {
    res <- t(t(x) * 1e6 / colSums(x))
  }
  else if(method == "RPKM") {
    if(length(len) != nrow(x)) stop('Length vector not matching gene number\n')
    res <- t(t(x) * 1e6 / colSums(x)) / (len / 1e3)
  }
  else if(method == "TPM") {
    if(length(len) != nrow(x)) stop('Length vector not matching gene number\n')
    x <- x / (len / 1e3)
    res <- t(t(x) * 1e6 / colSums(x))
  }
  else stop('This method is not available\n')
  
  return(as.data.frame(res))
}


#' Transforms data according to column-wise grouping using a function (e.g. mean, median...)
#' @param x A dataframe or matrix
#' @param group.lbls Vector with grouping of columns
#' @param FUN Function to be used for data transformation (e.g. function(x) apply(x,1,mean))
#' @return Data frame with transformed data
#' 
#' @author carlga
#'
groupTransform <- function(x, group.lbls, FUN) {
  
  group.lbls.uniq <- unique(group.lbls)
  group.lbls.uniq <- split(group.lbls.uniq, 1:length(group.lbls.uniq))
  
  res <- lapply(group.lbls.uniq, function(lbl) FUN(x[, group.lbls==lbl]))
  res <- dplyr::bind_cols(res)
  res <- as.data.frame(res)
  row.names(res) <- row.names(x)
  colnames(res) <- unlist(group.lbls.uniq)
  
  return(res)
}


#' Performs transcriptome-wide signal-to-noise ratio (tSNR)
#' @param x A dataframe or matrix
#' @param group.lbls Vector with grouping of columns
#' @return Data frame with tSNR values between groups
#' 
#' @author carlga
#'
tSNR <- function(x, group.lbls) {
  
  group.lbls.uniq <- unique(group.lbls)
  group.means <- groupTransform(x, 
                                group.lbls = group.lbls, 
                                FUN = function(x) apply(x,1,mean))
  
  comb <- expand.grid(group.lbls.uniq, group.lbls.uniq)
  comb <- split(comb, 1:nrow(comb))
  
  res <- lapply(comb, function(lbls) {
    
    # variance of X transcriptomes
    names.x <- names(x)[group.lbls==lbls[[1]]]
    names.x <- split(names.x, 1:length(names.x))
    dists.x <- lapply(names.x, function(i) {
      d <- dist(rbind(x[[i]], group.means[[lbls[[1]]]]), method = 'euclidean')
      return(d)
    })
    dists.x <- unlist(dists.x)
    var.x <- sum(dists.x^2)/(length(dists.x)-1)
    
    # variance of Y transcriptomes
    names.y <- names(x)[group.lbls==lbls[[2]]]
    names.y <- split(names.y, 1:length(names.y))
    dists.y <- lapply(names.y, function(j) {
      d <- dist(rbind(x[[j]], group.means[[lbls[[2]]]]), method = 'euclidean')
      return(d)
    })
    dists.y <- unlist(dists.y)
    var.y <- sum(dists.y^2)/(length(dists.y)-1)
    
    # calculate tSNR
    dist.x.y <- dist(rbind(group.means[[lbls[[1]]]], group.means[[lbls[[2]]]]), method = 'euclidean')
    tSNR <- dist.x.y/sqrt(var.x/length(dists.x)+var.y/length(dists.y))
    
    return(tSNR)
  })
  
  res <- unlist(res)
  res <- matrix(res, ncol = length(group.lbls.uniq))
  row.names(res) <- group.lbls.uniq
  colnames(res) <- group.lbls.uniq
  
  return(as.data.frame(res))
}


#' Function to load GMT files
#' @param file Path to file to be read
#' @return List of sets contained in the GMT file
#' 
#' @author Minghui Wang
#'
readGMT <- function(file) {
  if(!file.exists(file)) stop('File ',file,' not available\n')
  x <- readLines(file)
  n <- length(x)
  res <- list(genesets = vector(mode = "list", length = n),
              geneset.names = vector(mode = "character", length = n),
              geneset.descriptions = vector(mode = "character", length = n))
  for(i in 1:n) {
    s <- strsplit(x[i],'\t')[[1]]
    res$genesets[[i]] <- s[-c(1:2)]
    res$geneset.names[i] <- s[1]
    res$geneset.descriptions[i] <- s[2]
  }
  names(res$genesets) <- res$geneset.names
  res
}


#' Performs gene set enrichment analysis
#' @param gene.rnk A named vector of scores ranking genes
#' @param gene.sets List with gene sets to test enrichment
#' @param min.size Minimum size of gene sets to test
#' @param max.size Maximum size of gene sets to test
#' @param eps Sets boundary for calculating P value
#' @param nproc Sets BPPARAM to use nproc workers 
#' @return List with GSEA results
#' 
#' @author carlga
#'
runGSEA <- function(gene.rnk, gene.sets, min.size, max.size, eps = 0.0, nproc = 0) {
  require(fgsea)
  
  res <- list()
  res$all <- fgsea::fgsea(stats = gene.rnk,
                          pathways = gene.sets,
                          minSize = min.size,
                          maxSize = max.size,
                          eps = eps,
                          nproc = nproc)
  
  res$collapsed <- fgsea::collapsePathways(fgseaRes = res$all[order(pval)],
                                           pathways = gene.sets,
                                           stats = gene.rnk)
  res$collapsed <- res$all[pathway %in% res$collapsed$mainPathways]
  
  return(res)
}


#' Calculates similarity between two sets using selected method
#' @param x Character vector with elements from first set
#' @param y Character Vector with elements from second set
#' @param method Method used for determining similarity between sets
#' @return Similarity value
#' 
#' @author carlga
#'
getSimilarity <- function(x, y, method = NULL) {
  
  if(is.null(method)) stop('Please, indicate a valid method\n')
  
  set.x <- unique(x)
  set.y <- unique(y)
  
  if(method == "overlap") {
    intersection_n <- length(dplyr::intersect(set.x, set.y))
    min_n <- min(c(length(set.x), length(set.y)))
    res <- intersection_n/min_n
  }
  else if(method == "jaccard") {
    intersection_n <- length(intersect(set.x, set.y))
    union_n <- length(set.x) + length(set.y) - intersection_n
    res <- intersection_n/union_n
  }
  else stop('This method is not available\n')
  
  return(res)
}


#' Generates network from gene set enrichment analysis (GSEA) results
#' @param gsea.res List of GSEA results
#' @param gene.sets List with gene sets
#' @param gsea.cutoff Cutoff for determining significantly changed gene sets
#' @param similarity.method Approach for similarity calculation between gene sets 
#' ('overlap', 'jaccard' or 'combined')
#' @param similarity.cutoff Cutoff for filtering edges based on similarity
#' @return List with node and edge tables for similarity network of enriched gene sets
#' 
#' @author carlga
#'
buildGSEANetwork <- function(gsea.res, 
                             gene.sets, 
                             gsea.cutoff=0.05, 
                             similarity.method='combined', 
                             similarity.cutoff=0.5) {
  
  enriched <- lapply(gsea.res, function(x) x$all[padj<gsea.cutoff]$pathway)
  enriched <- unique(unlist(enriched))
  enriched <- split(enriched, 1:length(enriched))
  
  nodes <- lapply(enriched, function(set) {
    lapply(names(gsea.res), function(comp) {
      data.frame(id = rep(set, 3),
                 type = rep('ENR', 3),
                 genes = rep(paste(gene.sets[[set]], collapse = ','), 3),
                 size = rep(length(gene.sets[[set]]), 3),
                 name = paste(rep(comp, 3), c('NES', 'pval', 'padj'), sep = '_'),
                 value = unlist(gsea.res[[comp]]$all[pathway == set, c('NES', 'pval', 'padj')]))
    })
  })
  nodes <- dplyr::bind_rows(nodes)
  nodes <- tidyr::pivot_wider(nodes)
  
  edges <- lapply(enriched, function(set1) {
    edges <- lapply(enriched, function(set2) {
      if(similarity.method == 'combined') {
        overlap.sim <- getSimilarity(gene.sets[[set1]],gene.sets[[set2]], method = 'overlap')
        jaccard.sim <- getSimilarity(gene.sets[[set1]],gene.sets[[set2]], method = 'jaccard')
        sim <- (overlap.sim+jaccard.sim)/2
      }
      else {
        sim <- getSimilarity(gene.sets[[set1]],gene.sets[[set2]], method = similarity.method)
      }
      data.frame(source = set1,
                 target = set2,
                 interaction = 'overlap',
                 genes = paste(dplyr::intersect(gene.sets[[set1]], gene.sets[[set2]]), collapse = ','),
                 size = length(dplyr::intersect(gene.sets[[set1]], gene.sets[[set2]])),
                 similarity = sim)
    })
    dplyr::bind_rows(edges)
  })
  edges <- dplyr::bind_rows(edges)
  edges <- dplyr::filter(edges, source != target & similarity > similarity.cutoff)
  
  net <- list(nodes = nodes,
              edges = edges)
  
  return(net)
}


#' Calculates pathway activity scores (PAS) using pagoda2 method
#' @param x Seurat object with single cell data
#' @param gene.sets List with gene sets
#' @param npcs Number of principal components to use for dimensionality reduction
#' @return Data frame with single cell PAS values
#' 
#' @author carlga
#'
getPAS <- function(x, gene.sets, npcs = 50) {
  require(Seurat)
  require(pagoda2)
  
  mtx <- as.matrix(Seurat::GetAssayData(x, assay = "RNA", slot = "counts"))
  res <- x@meta.data
  
  cellnames <- colnames(mtx)
  genenames <- rownames(mtx)
  p2 <- Pagoda2$new(as.matrix(mtx), n.cores = 1, log.scale=F)
  
  p2$adjustVariance(plot=F)
  
  p2$calculatePcaReduction(nPcs = npcs, use.odgenes=F,fastpath=F)
  
  path_names <- c()
  env <- new.env(parent=globalenv())
  invisible(lapply(1:length(gene.sets),function(i) {
    genes <- intersect(gene.sets[[i]],genenames)
    name <- paste0(names(gene.sets[i]),i)
    if(length(genes)>3){
      assign(name, genes, envir = env)
      path_names <- c(path_names, name)
    }
  }))
  
  p2$testPathwayOverdispersion(setenv = env, verbose = T,
                               recalculate.pca = T,
                               min.pathway.size = 1, 
                               max.pathway.size = 2000)
  
  path_names <- gsub("^\\d+|\\d+$", "", names(p2$misc$pwpca))
  score <- matrix(NA,nrow=length(path_names),ncol=length(cellnames))
  rownames(score) <- path_names
  colnames(score) <- cellnames
  for(i in 1:length(p2$misc$pwpca)){
    if(!is.null(p2$misc$pwpca[[i]]$xp$score)){
      score[i,] <- as.numeric(p2$misc$pwpca[[i]]$xp$scores)
    }
  }
  score <- as.data.frame(t(score))
  res <- cbind(res, score)
  
  return(res)
}


#' Get percentage of distribution left or right to a given value
#' @param value A numeric value to check in the distribution
#' @param distr Vector of numeric values
#' @param tail String indicating which tail to calculate
#' @return Percentage of distribution left or right of value
#' 
#' @author carlga
#' 
tailFraction <- function(value, distr, tail = 'left') {
  if(tail == 'left') {
    res <- length(distr[distr <= value])/length(distr)
  }
  else if(tail == 'right'){
    res <- length(distr[distr >= value])/length(distr)
  }
  else {
    stop('Indicate tail as left or right')
  }
  
  return(res)
}


#' Makes duplicated elements in vector unique
#' @param v A vector with duplicated elements to make unique
#' @param sep Separator used to add unique number to duplicated elements
#' @return Uniquified vector
#' 
#' @author carlga
#' 
uniquify <- function(v, sep = '_'){
  ave(v, v, FUN = function(x) {
    if(length(x) > 1) {
      paste(x, 1:length(x), sep = sep)
    } 
    else x
  })
}
