#' Quantile Normalization Function
#'
#' @param x A numeric vector. The reference vector whose quantiles will be used for normalization.
#' @param y A numeric vector. The vector to be normalized based on the quantiles of `x`.
#'
#' @return A numeric vector of the same length as `y`, containing the quantile-normalized values.
#'
#' @importFrom stats approx
#' @export
zp.quantile <- function(x, y) {
  o.x <- order(x)
  o.y <- order(y)
  r.y <- rank(y, ties.method = "average")
  x <- x[o.x]
  y <- y[o.y]
  z.x <- seq(0, 1, length.out = length(x))
  z.y <- seq(0, 1, length.out = length(y))
  new.y <- stats::approx(x = z.x, y = x, xout = z.y)$y
  y <- new.y[r.y]
  return(y)
}

#' Move Outliers Function
#'
#' @param x A numeric vector or matrix.
#' @param d The distance threshold for outlier detection. If NULL, it's calculated from the data.
#' @param fraction The fraction of points to consider as outliers.
#'
#' @return A numeric vector or matrix with outliers moved.
#'
#' @importFrom DDoutlier DB
#' @importFrom stats quantile
move.outlier <- function(x, d = NULL, fraction = 0.01) {
  x.is.vector <- is.null(dim(x))
  if (x.is.vector) x <- matrix(x, ncol = 1)
  
  for (n in 1:ncol(x)) {
    this.x <- x[, n]
    if (is.null(d)) {
      this.d <- diff(quantile(this.x, c(0.25, 0.75))) / 3
    } else {
      this.d <- d
    }
    
    this.out <- DDoutlier::DB(matrix(this.x, ncol = 1), this.d, fraction)
    this.sel <- which(this.out$classification == "Outlier")
    
    if (length(this.sel) > 0) {
      new.x <- zp.quantile(this.x[-this.sel], this.x)
      x[, n] <- new.x
    }
  }
  
  if (x.is.vector) x <- as.vector(x)
  return(x)
}

#' Multidimensional Embedding Visualization (ADM)
#'
#' @param e A list of matrices, each representing a dimension reduction result.
#' @param k.dim Number of dimensions to use. If NULL, uses the number of columns in the first matrix of e.
#' @param dist.power Power to raise distances to when merging. Default is 0.5.
#' @param conn.prop Proportion of connections to consider. Default is 0.02.
#' @param raw.d.pwr Power for raw distance calculation. Default is 0.5.
#' @param diffu.steps Steps for diffusion. If NA, calculated based on graph properties.
#' @param diffu.factor Factor to multiply diffusion steps by. Default is 3.5.
#' @param distr.template Distribution template to use. Options are "none", "combine", "gamma", or "parametric". Default is "gamma".
#' @param gamma.shape Shape parameter for gamma distribution. Default is 3.
#' @param gamma.rate Rate parameter for gamma distribution. Default is 3.
#' @param scale.dist Whether to scale distances. Default is TRUE.
#' @param symmetrize Method to symmetrize distance matrix. Options are "mean", "min", or "max". Default is "mean".
#' @param dist.quantile Quantile to use for distance calculations. Default is 0.25.
#'
#' @return A list containing the final merged distance matrix.
#'
#' @importFrom stats dist quantile cor rgamma
#' @importFrom utils write.csv
#' @importFrom diffudist get_distance_matrix_from_T
#' @importFrom igraph graph_from_adjacency_matrix shortest.paths
#' @importFrom BiocParallel bplapply
#' @importFrom fitdistrplus fitdist
#'
#' @export
adm <- function(e, k.dim = NULL, dist.power = 0.5, conn.prop = 0.02, raw.d.pwr = 0.5, 
                diffu.steps = NA, diffu.factor = 3.5, distr.template = "gamma", 
                gamma.shape = 3, gamma.rate = 3, scale.dist = TRUE, 
                symmetrize = "mean", dist.quantile = 0.25) {
  
  if (is.null(k.dim)) k.dim <- ncol(e[[1]])
  print("working on ADM...")
  fake.fun <- function(this.e, raw.d.pwr, diffu.steps, conn.prop, scale.dist, symmetrize, diffu.factor, dist.quantile) {
    this.e <- move.outlier(this.e)
    this.d <- as.matrix(stats::dist(this.e))
    n <- nrow(this.d)
    
    diag(this.d) <- Inf
    
    conn.cutoff <- stats::quantile(as.dist(this.d), conn.prop)
    this.conn <- 1 * (this.d <= conn.cutoff)
    
    for (i in 1:nrow(this.d)) {
      sel <- which(this.d[i,] <= stats::quantile(this.d[i,], conn.prop))
      this.conn[i, sel] <- 1
    }
    
    this.graph <- igraph::graph_from_adjacency_matrix(this.conn, mode = "undirected")
    
    pmat <- 1 / this.d^raw.d.pwr
    diag(pmat) <- 0
    pmat[pmat == Inf] <- max(pmat[pmat != Inf], na.rm = TRUE)
    pmat[is.na(pmat)] <- max(pmat, na.rm = TRUE)
    
    pmat <- pmat * this.conn
    for (i in 1:nrow(pmat)) pmat[i,] <- pmat[i,] / sum(pmat[i,])
    
    if (is.na(diffu.steps[1])) {
      sp.mat <- igraph::shortest.paths(this.graph)
      sp.mat[sp.mat == Inf] <- NA
      mean.steps <- apply(sp.mat, 1, stats::quantile, na.rm = TRUE, probs = dist.quantile)
      diffu.steps <- stats::quantile(mean.steps, c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE) * diffu.factor
    }
    
    all.d <- vector("list", length(diffu.steps))
    for (i in seq_along(diffu.steps)) {
      this.d <- as.matrix(suppressMessages(diffudist::get_distance_matrix_from_T(pmat, diffu.steps[i])))
      if (scale.dist) {
        this.d <- 1 - stats::cor(this.d, method = "spearman")
        this.d[this.d > 1] <- 1
      }
      all.d[[i]] <- this.d
    }
    
    merge.d <- all.d[[1]]
    for (i in 1:nrow(merge.d)) {
      this.dist <- abs(mean.steps[i] * diffu.factor - diffu.steps)
      this.closest <- which.min(this.dist)[1]
      merge.d[i,] <- all.d[[this.closest]][i,]
    }
    this.d <- (merge.d + t(merge.d)) / 2
    
    if (symmetrize == "mean") {
      this.d <- this.d + t(this.d)
    } else if (symmetrize == "min") {
      this.d.2 <- t(this.d)
      this.d[this.d > this.d.2] <- this.d.2[this.d > this.d.2]
    } else if (symmetrize == "max") {
      this.d.2 <- t(this.d)
      this.d[this.d < this.d.2] <- this.d.2[this.d < this.d.2]
    }
    
    list(this.d = this.d)
  }
  
  d <- BiocParallel::bplapply(e, fake.fun, raw.d.pwr = raw.d.pwr, diffu.steps = diffu.steps, 
                              conn.prop = conn.prop, scale.dist = scale.dist, 
                              symmetrize = symmetrize, diffu.factor = diffu.factor, 
                              dist.quantile = dist.quantile)
  
  if (distr.template != "none") {
    if (distr.template == "combine") {
      d.template <- unlist(lapply(d, function(x) as.dist(x[[1]])))
    } else if (distr.template == "gamma") {
      all.fit <- t(sapply(d, function(x) {
        fitdistrplus::fitdist(as.vector(as.dist(x[[1]])), "gamma", method = "mme")$estimate
      }))
      ave.fit <- colMeans(all.fit)
      d.template <- stats::rgamma(nrow(d[[1]][[1]]) * ncol(d[[1]][[1]]), shape = ave.fit[1], rate = ave.fit[2])
    } else if (distr.template == "parametric") {
      d.template <- stats::rgamma(nrow(d[[1]][[1]]) * ncol(d[[1]][[1]]), shape = gamma.shape, rate = gamma.rate)
    }
    
    for (m in seq_along(d)) {
      this.d <- as.dist(d[[m]][[1]])
      this.d2 <- zp.quantile(d.template, this.d)
      attributes(this.d2) <- attributes(this.d)
      this.d2 <- as.matrix(this.d2)
      d[[m]][[1]] <- this.d2
    }
  }
  
  dd <- Reduce(`+`, lapply(d, function(x) x[[1]]^dist.power))
  
  list(diffu.dist = dd)
}


