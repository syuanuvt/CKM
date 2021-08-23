#' The function extends the existing function for calculating the Gap statistics (i.e., clusGap in the package Cluster) such that users can also specify the minimal value of K (i.e. the number of clusters)
#' The other paramters as well as the functions are identical to the original function
#' See ?clusGap for more information
#'
#' @param dataset  a numeric matrix or a data frame
#' @param FUNcluster the main function of clustering analysis (see ?clusGap for more information)
#' @param minclust the minimal possible number of clusters
#' @param maxclust the maximal possible number of clusters
#' @param B        the number of bootstrap samples
#' @param d.power  a positive integer specifying the power that is appled to the euclidean distances
#' @param spaceH0  either "scaledPCA" or "original" that specifies the space of the H0 distribution
#' @param verbose  a logical value determining whether the process of computation is shown
#' @param nstart.o number of resamples for the original sample
#' @param nstart.b number of resamples for the bootstrapping sample
#' @return The function will return a ckm object that is the list of five elements. The first denotes the selected number of masking variables; the second includes all indicies of signaling variables; the third is a vector illustrating cluster assignment; the forth is the pre-determined or selected "optimal" number of clusters; the fifth is the original dataset.
#' @examples
#' maxcluster <- 10
#' ckm.sel.all <- CKMSelAll(dataset, maxcluster)
#'
clusGap_con <- function (x, FUNcluster, K.min, K.max, B = 100, d.power = 1, spaceH0 = c("scaledPCA", "original"), verbose = interactive(), nstart.o = 25, nstart.b = 25)
{
  stopifnot(is.function(FUNcluster), length(dim(x)) == 2,
            K.max >= 2, (n <- nrow(x)) >= 1, ncol(x) >= 1)
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0)
    stop("'B' has to be a positive integer")
  cl. <- match.call()
  if (is.data.frame(x))
    x <- as.matrix(x)
  ii <- seq_len(n)
  W.k <- function(X, kk, nstart.o) {
    clus <- if (kk > 1)
      FUNcluster(X, kk, nstart = nstart.o)$cluster
    else rep.int(1L, nrow(X))
    0.5 * sum(vapply(split(ii, clus), function(I) {
      xs <- X[I, , drop = FALSE]
      sum(dist(xs)^d.power/nrow(xs))
    }, 0))
  }
  W.k.b <- function(X, kk, nstart.b) {
    clus <- if (kk > 1)
      FUNcluster(X, kk, nstart.b)$cluster
    else rep.int(1L, nrow(X))
    0.5 * sum(vapply(split(ii, clus), function(I) {
      xs <- X[I, , drop = FALSE]
      sum(dist(xs)^d.power/nrow(xs))
    }, 0))
  }
  searchClust <- length(K.min:K.max)
  logW <- E.logW <- SE.sim <- numeric(searchClust)
  if (verbose)
    cat("Clustering k = K.min", K.min,",..., K.max (= ", K.max, "): .. ",
        sep = "")
  for (k in K.min:K.max) logW[k-K.min+1] <- log(W.k(x, k, nstart.o))
  if (verbose)
    cat("done\n")
  spaceH0 <- match.arg(spaceH0)
  xs <- scale(x, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  switch(spaceH0, scaledPCA = {
    V.sx <- svd(xs, nu = 0)$v
    xs <- xs %*% V.sx
  }, original = {
  }, stop("invalid 'spaceH0':", spaceH0))
  rng.x1 <- apply(xs, 2L, range)
  logWks <- matrix(0, B, searchClust)
  if (verbose)
    cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n",
        sep = "")
  for (b in 1:B) {
    z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1],
                                                 max = M[2]), nn = n)
    z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx),
                original = z1) + m.x
    for (k in K.min:K.max) {
      logWks[b, (k-K.min+1)] <- log(W.k.b(z, k, nstart.b))
    }
    if (verbose)
      cat(".", if (b%%50 == 0)
        paste(b, "\n"))
  }
  if (verbose && (B%%50 != 0))
    cat("", B, "\n")
  E.logW <- colMeans(logWks)
  SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
  structure(class = "clusGap", list(Tab = cbind(logW, E.logW,
                                                gap = E.logW - logW, SE.sim), call = cl., spaceH0 = spaceH0,
                                    n = n, B = B, FUNcluster = FUNcluster))
}
