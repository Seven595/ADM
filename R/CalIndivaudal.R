#' @importFrom stats dist quantile cmdscale
#' @importFrom utils globalVariables
#' @importFrom rARPACK svds
#' @importFrom MASS isoMDS sammon
#' @importFrom uwot umap
#' @importFrom phateR phate
NULL

utils::globalVariables(c("Rtsne"))

#' Perform Multiple Dimensionality Reduction Methods
#'
#' @description
#' This function applies various dimensionality reduction methods to the input data.
#' It supports methods such as PCA, MDS, iMDS, Sammon mapping, LLE, HLLE, Isomap,
#' kPCA, LEIM, UMAP, t-SNE, PHATE, and KEF.
#'
#' @param data A numeric matrix or data frame where rows are observations and columns are features.
#' @param dim Integer. The number of dimensions for the reduced space. Default is 2.
#' @param methods Character vector. Methods to be applied. Default includes all supported methods.
#' @param kpca.sigma Numeric vector. Sigma values for kPCA. Default is c(0.001, 0.002).
#' @param umap.k Integer vector. Number of neighbors for UMAP. Default is c(30, 50).
#' @param tsne.perplexity Numeric vector. Perplexity values for t-SNE. Default is c(30, 50).
#' @param phate.k Integer vector. Number of nearest neighbors for PHATE. Default is c(30, 50).
#' @param cal_dist Logical. Whether to calculate distance matrix. Default is TRUE.
#'
#' @return A list containing two elements:
#'   \item{embed.list}{A list of matrices, each representing the reduced data for a method}
#'   \item{method_name}{A character vector of method names corresponding to embed.list}
#'
#' @details
#' The function applies each specified method to the input data. For some methods
#' (kPCA, UMAP, t-SNE, PHATE), multiple parameter settings are tried, resulting in
#' multiple outputs for these methods.
#'
#' @note
#' - This function requires several packages including rARPACK, MASS, dimRed, umap, and phateR.
#' - The function prints progress messages to the console.
#' - Some methods may fail if the required packages are not installed.
#'
#' @examples
#' \dontrun{
#' # Assuming 'data' is your input dataset
#' result <- candidate.visual(data, dim = 2,
#'                            methods = c("PCA", "MDS", "UMAP"),
#'                            umap.k = c(15, 30))
#' pca_result <- result$embed.list[[1]]
#' method_names <- result$method_name
#' }
#'
#' @export
candidate.visual <- function(data, dim=2, methods= c("PCA", "MDS", "iMDS", "Sammon", "HLLE","Isomap",
                                                     "kPCA", "LEIM", "UMAP", "tSNE", "PHATE","KEF"),
                             kpca.sigma = c(0.001, 0.002),
                             umap.k= c(30, 50),
                             tsne.perplexity = c(30, 50),
                             phate.k = c(30,50),
                             cal_dist = TRUE){

  n=dim(data)[1]
  dim.red.data = list()
  if(cal_dist){
    dist.data = stats::dist(data)
  }
  print(dim(data))
  name.method =c()
  ####################
  ###### PCA
  ####################

  i=0
  if(sum(methods == "PCA")>0){
    i=i+1
    print("PCA calculating...")
    pc.data = rARPACK::svds(as.matrix(data), k =dim)
    dim.red.data[[i]] = pc.data$u[,1:2]
    name.method =  c(name.method, "PCA")
  }


  #####################
  ####### classical MDS
  #####################

  if(sum(methods == "MDS")>0){
    i=i+1
    print("MDS calculating...")
    dim.red.data[[i]] = stats::cmdscale(dist.data, k=dim)
    name.method =  c(name.method, "MDS")
  }
  #####################
  ####### isoMDS
  #####################

  if(sum(methods == "iMDS")>0){
    i=i+1
    print("iMDS calculating...")
    imds.data = MASS::isoMDS(dist.data, k=dim)
    dim.red.data[[i]] = imds.data$points
    name.method =  c(name.method, "iMDS")
  }

  #####################
  ####### Sammon's nonlinear mapping
  #####################

  if(sum(methods == "Sammon")>0){
    i=i+1
    print("Sammon calculating...")
    sam.data = MASS::sammon(dist.data, k=dim)
    dim.red.data[[i]] = sam.data$points
    name.method =  c(name.method, "Sammon")
  }

  #####################
  ####### HLLE
  #####################

  if(sum(methods == "HLLE")>0){
    i=i+1
    print("HLLE calculating...")
    if(requireNamespace("dimRed", quietly = TRUE)) {
      hlle.data <- dimRed::embed(data, "HLLE", knn = 20, ndim = dim)
      dim.red.data[[i]] = hlle.data@data@data
      name.method =  c(name.method, "HLLE")
    } else {
      warning("dimRed package is not available. Skipping HLLE calculation.")
    }
  }

  #####################
  ####### isomap
  #####################

  if(sum(methods == "Isomap")>0){
    i=i+1
    print("Isomap calculating...")
    if(requireNamespace("dimRed", quietly = TRUE)) {
      imp.data <- dimRed::embed(data, "Isomap", knn = 20, ndim = dim)
      dim.red.data[[i]] = imp.data@data@data
      name.method =  c(name.method, "Isomap")
    } else {
      warning("dimRed package is not available. Skipping Isomap calculation.")
    }
  }

  #####################
  ####### kPCA
  #####################

  if(sum(methods == "kPCA")>0){
    for(j in 1:length(kpca.sigma)){
      i=i+1
      print("kPCA calculating...")
      if(requireNamespace("dimRed", quietly = TRUE)) {
        kpca.data <- dimRed::embed(data, "kPCA", kpar = list(sigma = kpca.sigma[j]), ndim = dim)
        dim.red.data[[i]] = kpca.data@data@data
        name.method =  c(name.method, paste0("kPCA", j))
      } else {
        warning("dimRed package is not available. Skipping kPCA calculation.")
      }
    }
  }

  #####################
  ####### Laplacian Eigenmap
  #####################

  if(sum(methods == "LEIM")>0){
    i=i+1
    print("LEIM calculating...")
    if(requireNamespace("dimRed", quietly = TRUE)) {
      lem.data <- dimRed::embed(data, "LaplacianEigenmaps", ndim = dim)
      dim.red.data[[i]] = lem.data@data@data
      name.method =  c(name.method, "LEIM")
    } else {
      warning("dimRed package is not available. Skipping LEIM calculation.")
    }
  }

  #####################
  ####### UMAP
  #####################

  if(sum(methods == "UMAP")>0){
    for(j in 1:length(umap.k)){
      i=i+1
      print("UMAP calculating...")
      umap.data <- uwot::umap(data,  n_neighbors = umap.k[j], n_components = dim)
      dim.red.data[[i]] = umap.data
      name.method =  c(name.method, paste0("UMAP", j))
    }
  }
  #####################
  ####### tSNE
  #####################

  if(sum(methods == "tSNE")>0){
    for(j in 1:length(tsne.perplexity)){
      i=i+1
      print("tSNE calculating...")
      if(requireNamespace("Rtsne", quietly = TRUE)) {
        tsne.data <- Rtsne::Rtsne(data, perplexity = tsne.perplexity[j], dims = dim)
        dim.red.data[[i]] = tsne.data$Y
        name.method =  c(name.method, paste0("tSNE", j))
      } else {
        warning("Rtsne package is not available. Skipping t-SNE calculation.")
      }
    }
  }

  #####################
  ####### PHATE
  #####################

  if(sum(methods == "PHATE")>0){
    for(j in 1:length(phate.k)){
      i=i+1
      print("PHATE calculating...")
      dim.red.data[[i]]  <- phateR::phate(data, knn = phate.k[j], ndim = dim)$embedding
      name.method =  c(name.method, paste0("PHATE", j))
    }
  }

  #####################
  ####### KEF
  #####################
  if(sum(methods == "KEF") > 0){
    i = i + 1
    print("KEF calculating...")
    dist.mat = stats::dist(data)
    K.mat = exp(-as.matrix(dist.mat)^2/stats::quantile(dist.mat,0.5)^2)
    eigen.K.mat = eigen(K.mat)
    kef_u1 = eigen.K.mat$vectors[,2] * eigen.K.mat$values[2]
    kef_u2 = eigen.K.mat$vectors[,3] * eigen.K.mat$values[3]
    print(dim(kef_u1))
    dim.red.data[[i]] = cbind(kef_u1, kef_u2)
    name.method =  c(name.method, "KEF")
  }

  return(list(embed.list=dim.red.data, method_name=name.method))
}