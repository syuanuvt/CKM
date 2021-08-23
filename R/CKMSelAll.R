#' The function is an advanced version of CKMSelVar that selects the number of clusters together with the number of masking variables
#'
#' @param dataset the orginal dataset on which CKM and its model selection procedure operates
#' @param minclust the minimal possible number of clusters
#' @param maxclust the maximal possible number of clusters
#' @param search the mode of selecting over the grid. "all" = selecting over each point of the grid; while it maximizes the accuracy, it is overly slow with large number of variables.
#' "sub" = the "grid search with a zoom" strategy; while it is less accurate compared to searching the full grid, it is efficient even with large number of variables.
#' "dep" automatically adjust to one of the above two methods based on the number of variables. When # variables < 25, the search covers every possible value of the grid. This is also the default option.
#' @param maxnum the parameter is only useful when the "grid search with a zoom" strategy is applied. It restricts the maximal number of values searched over in any iteration. The default value is set at 10.
#' @param n.rep the number of permutated datasets when calculating the gap statistic
#' @param method different criterion exists as to determine the number of clusters based on the gap statistic; the users can try out these various options in the "method" argument. The default is "globalmax": selects the largest gap over all possible number of clusters (i.e. global maxima). Other options include "firstSEmax": select the "first" gap that falls within the range of the largest gap minus one SE (i.e. the one SE role); "firstmax": select the first largest gap (i.e., local maxima), "Tibs2001SEmax": the recommened guideline of Tibshirani, 2011 that takes the one-SE rule.
#' @param kmeans_starts the number of starts used in the kmeans algorithm
#' @param recal whether a final step of KM is carried out to update cluster partitions (recommended when the number of signaling variable is small)
#' @param sr whether a scree ratio test has been added in the end (recommended when the number of clusters k is large (i.e., k > 15))
#' @return The function will return a ckm object that is the list of five elements. The first denotes the selected number of masking variables; the second includes all indicies of signaling variables; the third is a vector illustrating cluster assignment; the forth is the pre-determined or selected "optimal" number of clusters; the fifth is the original dataset.
#' @examples
#' mincluster <- 2
#' maxcluster <- 6
#' ckm.sel.all <- CKMSelAll(dataset, minclust, maxcluster)
#'

CKMSelAll <- function(dataset, minclust, maxclust, search = "dep", maxnum = 10, n.rep = 20, method = "globalmax", kmeans_starts = 20, recal = TRUE, sr = TRUE){

  n.obs <- nrow(dataset)
  n.var <- ncol(dataset)

  if(search == "dep"){
    search <- "sub"
    if(n.var < 25){
      search <- "all"
    }
  }

  results.sel <- list()
  partition.sel <- list()
  search.grid <- minclust:maxclust
  n.sel <- rep(NA, length(search.grid))

  for (i in minclust:maxclust){
    a <- try(CKMSelVar(dataset, n.cluster = i, search = search, maxnum = maxnum, n.rep = n.rep, kmeans_starts = kmeans_starts, recal, sr), silent = TRUE)
    if(!inherits(a,'try-error')){
      results.sel[[(i-minclust+1)]] <- a$signaling.set
      partition.sel[[(i-minclust+1)]] <- a$cluster.assign
      n.sel[(i-minclust+1)] <- a$n.noisevar
    }
    if(inherits(a,'try-error')){
      break
    }
  }

  if(i != minclust){
    if(inherits(a,'try-error')){
      maxclust <- (i-1)
    }
    valid.set.all <- Reduce(intersect, results.sel)
    new.data <- dataset[,valid.set.all]
    b <- clusGap_con(new.data, FUN=kmeans, K.min = minclust, K.max=maxclust, B=50, spaceH0 = "original",
                   nstart.o = 500, nstart.b = 30)
  #b <- clusGap_con(new.data, FUN=kmeans, nstart = n.rep, K.min = minclust, K.max=maxclust, B=boot.rep, spaceH0 = "original")
    opt.cluster <- maxSE(b$Tab[,3], b$Tab[,4], method = method)

    n.noisevar <- n.sel[opt.cluster]
    signaling.set <- results.sel[[opt.cluster]]
    cluster.assign <- partition.sel[[opt.cluster]]
    opt.cluster <- opt.cluster + minclust-1
    results <- list(n.noisevar = n.noisevar, signaling.set = signaling.set, cluster.assign = cluster.assign, opt.cluster = opt.cluster, org.data = dataset)
    class(results) <- "ckm"
  }
  if(i == minclust){
    results <- NULL
  }
  return(results)
}
