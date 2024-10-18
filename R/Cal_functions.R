#' Perform UMAP and calculate distance-based group statistics
#'
#' This function applies UMAP to the input data and calculates distance-based group statistics.
#'
#' @param x A distance matrix or object that can be coerced to one.
#' @param info A vector of group labels for each data point.
#' @param do.plot Logical; if TRUE, creates a pairs plot of the UMAP result.
#' @param n_components Integer; the number of dimensions for UMAP.
#' @param k Vector of integers; the number of neighbors to consider for distance calculations.
#'
#' @return A vector of mean distance-based group statistics across multiple UMAP runs.
#'
#' @importFrom uwot umap
#' @importFrom stats dist quantile
#' @importFrom graphics pairs
#' @importFrom grDevices rainbow
#' @export
umap5 <- function(x, info, do.plot=TRUE, n_components, k=c(1,2,5,10,20))
{
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is needed for this function to work. Please install it.", call. = FALSE)
  }

  rec <- NULL
  for(n in 1:5)
  {
    ensemble.data <- uwot::umap(stats::dist(x), n_components = n_components)
    ds <- dist.grp(as.matrix(stats::dist(ensemble.data)), info, k=k)[,2]
    if(n==1) rec <- ds
    else rec <- cbind(rec,ds)
  }
  if(do.plot) graphics::pairs(ensemble.data, col=grDevices::rainbow(length(unique(info)))[as.numeric(info)], pch=(1:length(unique(info)))[as.numeric(info)])
  to.return <- apply(rec,1,mean)
  return(to.return)
}


#' Calculate distance-based group statistics
#'
#' This function calculates distance-based group statistics for a given distance matrix and grouping.
#'
#' @param distmat A distance matrix.
#' @param grp A vector of group labels for each data point.
#' @param k Vector of integers; the number of neighbors to consider.
#'
#' @return A matrix with two columns: the k values and corresponding distance-based group statistics.
#'
#' @importFrom BiocParallel bplapply
#' @export
dist.grp <- function(distmat, grp, k=1:20)
{
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("Package 'BiocParallel' is needed for this function to work. Please install it.", call. = FALSE)
  }

  grpmat <- matrix(0,nrow=length(grp),ncol=length(grp))
  for(i in 1:length(grp)) for(j in 1:length(grp)) if(grp[i]==grp[j]) grpmat[i,j] <- 1

  rec <- cbind(k,k)
  diag(distmat) <- Inf

  dist.fun <- function(distmat, this.k, grpmat)
  {
    r <- NULL
    for(i in 1:nrow(distmat))
    {
      sel <- which(distmat[i,] <= stats::quantile(distmat[i,], this.k/ncol(distmat)))
      r <- c(r,grpmat[i,sel])
    }
    r <- sum(r==1)/length(r)
    r
  }

  rec[,2] <- unlist(BiocParallel::bplapply(k, dist.fun, distmat=distmat, grpmat=grpmat))
  rec
}

#' Calculate Cluster Conservation Index (CCI)
#'
#' This function calculates the Cluster Conservation Index (CCI) for both the ensemble and ADM methods.
#'
#' @param ensemble.out A list containing the output from the ensemble method, including the ensemble distance matrix.
#' @param mev.out A list containing the output from the ADM method, including the diffusion distance matrix.
#' @param info A vector or factor containing cluster information for each data point.
#'
#' @return A list containing two elements:
#'   \item{cci_viz}{A matrix with CCI values for the ensemble method}
#'   \item{cci_adm}{A matrix with CCI values for the ADM method}
#'
#' @details
#' The function calculates CCI for different numbers of neighbors (1, 2, 5, 10, 20) and
#' for different representations of the data (raw, UMAP, PCA).
#'
#' @importFrom stats dist eigen
#' @export
cal_cci <- function(ensemble.out, mev.out, info){
    n_components = 3
    cci_list = list()
    cci_viz <- dist.grp(ensemble.out$ensemble.dist.mat, info, k=c(1,2,5,10,20))
    cci_viz <- cbind(cci_viz, umap5(x=as.matrix(ensemble.out$ensemble.dist.mat), info=info,  n_components = n_components, do.plot=FALSE))
    cci_viz <- cbind(cci_viz, dist.grp(as.matrix(stats::dist(eigen(max(ensemble.out$ensemble.dist.mat)-ensemble.out$ensemble.dist.mat)$vec[,1:n_components])), info, k=c(1,2,5,10,20))[,2])
    colnames(cci_viz) <- c("n_neighbors","viz_raw","viz_umap","viz_pca")
    cci_list[[1]] <- cci_viz
    cci_adm <- dist.grp(mev.out$diffu.dist, info, k=c(1,2,5,10,20))
    cci_adm <- cbind(cci_adm, umap5(x=mev.out$diffu.dist, info=info,  n_components = n_components, do.plot=FALSE))
    cci_adm <- cbind(cci_adm, dist.grp(as.matrix(stats::dist(eigen(max(mev.out$diffu.dist)-mev.out$diffu.dist)$vec[,1:n_components])), info, k=c(1,2,5,10,20))[,2])
    colnames(cci_adm) <- c("n_neighbors","adm_raw","adm_umap","adm_pca")
    cci_list[[2]] <- cci_adm
    return (cci_list)
}

#' Calculate Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI)
#'
#' @param lowDim_data A matrix or data frame of low-dimensional data.
#' @param k The number of clusters for k-means clustering.
#' @param method_name A string naming the method being evaluated.
#' @param seed An integer seed for reproducibility.
#' @param info A vector of true class labels. If NULL, it will use the global 'info' variable.
#' @importFrom mclust adjustedRandIndex
#' @importFrom aricode NMI
#' @importFrom stats kmeans
#' @return A data frame with ARI and NMI scores.
#'
#' @details
#' This function performs k-means clustering on the input data and calculates
#' the Adjusted Rand Index and Normalized Mutual Information between the
#' clustering result and a pre-existing classification.
#'
#' @examples
#' data <- matrix(rnorm(200), 100, 2)
#' info <- sample(1:3, 100, replace = TRUE)  # Simulated cluster information
#' result <- cal_ari_nmi(data, k = 3, method_name = "Example Method", seed = 123, info = info)
#'
#' @export
cal_ari_nmi <- function(lowDim_data, k, method_name, seed, info = NULL){
  if (is.null(info)) {
    if (!exists("info")) {
      stop("'info' not provided and not found in global environment")
    }
    info <- get("info", envir = .GlobalEnv)
  }
  # set.seed(seed)
  cluster_viz <- stats::kmeans(lowDim_data, centers = k)
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("aricode", quietly = TRUE)) {
    stop("Package 'aricode' is needed for this function to work. Please install it.", call. = FALSE)
  }
  ARI <- mclust::adjustedRandIndex(info, cluster_viz$cluster)
  NMI <- aricode::NMI(info, cluster_viz$cluster)
  rec_ari_nmi <- data.frame(ARI = ARI, NMI = NMI)
  cat("******", method_name, "******\n")
  print(rec_ari_nmi)
  return(rec_ari_nmi)
}