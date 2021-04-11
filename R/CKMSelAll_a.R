###########################################################
##only used for the simulation, not included in the package
###########################################################

CKMSelAll_a <- function(dataset, maxclust, search = "dep", maxnum = 10, n.rep = 20, method = "globalmax", kmeans_starts = 10){

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
  n.sel <- rep(NA, (maxclust - 1))
  gap.sel <- rep(NA, (maxclust - 1))

  for (i in 2:maxclust){
    a <- CKMSelVar(dataset, n.cluster = i, search = search, maxnum = maxnum, n.rep = n.rep, kmeans_starts = kmeans_starts)
    results.sel[[(i-1)]] <- a$signaling.set
    partition.sel[[(i-1)]] <- a$cluster.assign
    n.sel[(i-1)] <- a$n.noisevar
    gap.sel[(i-1)] <- a$gap
  }

  valid.set.all <- Reduce(intersect, results.sel)
  new.data <- dataset[,valid.set.all]
  b <- clusGap(new.data, FUN=kmeans, K.max=maxclust, B=500, spaceH0 = "original")
  opt.cluster <- maxSE(b$Tab[,3], b$Tab[,4], method = method)

  n.noisevar <- n.sel[(opt.cluster - 1)]
  signaling.set <- results.sel[[(opt.cluster - 1)]]
  cluster.assign <- partition.sel[[(opt.cluster - 1)]]

  results <- list(stable.set = valid.set.all, cluster.assign = cluster.assign, opt.cluster = opt.cluster, b = b, gap = gap.sel, results.sel = results.sel)
  return(results)
}
