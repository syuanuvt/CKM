#' The function, combined with CardKMeans, selects the number of masking variables, given the fixed number of clusters
#'
#' @param dataset the orginal dataset on which CKM and its model selection procedure operates
#' @param n.cluster the total number of clusters
#' @param search the mode of selecting over the grid. "all" = selecting over each point of the grid; while it maximizes the accuracy, it is overly slow with large number of variables.
#' "sub" = the "grid search with a zoom" strategy; while it is less accurate compared to searching the full grid, it is efficient even with large number of variables.
#' "dep" automatically adjust to one of the above two methods based on the number of variables. When # variables < 25, the search covers every possible value of the grid. This is also the default option.
#' @param maxnum the parameter is only useful when the "grid search with a zoom" strategy is applied. It restricts the maximal number of values searched over in any iteration. The default value is set at 10.
#' @param n.rep the number of permutated datasets when calculating the gap statistic
#' @param kmeans_starts the number of starts used in the kmeans algorithm
#' @return @return The function will return a ckm object that is the list of five elements. The first denotes the selected number of masking variables; the second includes all indicies of signaling variables; the third is a vector illustrating cluster assignment; the forth is the pre-determined or selected "optimal" number of clusters; the fifth is the original dataset.
#' @examples
#' ncluster <- 3
#' ckm.sel.var <- CKMSelVar(dataset, ncluster)
#'

CKMSelVar <- function(dataset, n.cluster, search = "dep", maxnum = 10, n.rep = 20, kmeans_starts = 10){

  n.obs <- nrow(dataset)
  n.var <- ncol(dataset)

  if(search == "dep"){
    search <- "sub"
    if(n.var < 25){
      search <- "all"
    }
  }
  ## create the permutations: n.rep number of datasets
  permutate.data <- list()
  for (i in 1:n.rep){
    permutate.data[[i]] <- matrix(nrow = n.obs, ncol = n.var)
    for (j in 1:n.var){
      permutate.data[[i]][,j] <- dataset[,j][base::sample(1:n.obs, size = n.obs, replace = FALSE)]
    }
  }

  if(search == "all"){
    search.grid <- 1:(n.var-2-floor(n.var / 500))
    gap <-  rep(NULL, n.var)

    ## the original dataset
    for (i in 1:(n.var-2)) {
      method3.card <- CardKMeans(dataset, n.cluster, search.grid[i], kmeans_starts)
      valid.set <- method3.card$variables
      ori.between <- method3.card$betss

      rep.bet <- rep(NULL, n.rep)
      for (j in 1:n.rep){
        delta.data <- permutate.data[[j]][,valid.set]
        rep.bet[j] <- kmeans(delta.data, n.cluster, nstart = 5, iter.max = 50)$betweenss
      }

      gap[i] <- log(ori.between) - mean(log(rep.bet))
    }
    opt.n.noise <- which.max(gap)
  }

  if(search == "sub"){

    stop <- 0
    start.pos <- 1
    final.pos <- (n.var - 2)
    iter <- 0
    grid.test <- list()
    gap.all <- list()

    while(stop == 0){

      iter <- iter + 1
      step <- (final.pos - start.pos) / (maxnum - 1)

      ## the last round (step = 1)
      if (step  <= 1){
        step <- 1
        stop <- 1
        grid.iter <- seq(from = start.pos, by = step, to = final.pos)

      }
      ## not the last round (step != 1)
      else {
        step <- floor(step)
        grid.iter <- floor(seq(from = start.pos, to = final.pos, length.out = maxnum))
      }

      grid.test[[iter]] <- grid.iter

      gap <- rep(NULL, length(grid.iter))
      for (i in 1:length(grid.iter)){
        method3.card <- CardKMeans(dataset, n.cluster, grid.iter[i], kmeans_starts)
        valid.set <- method3.card$variables
        ori.between <- method3.card$betss
        rep.bet <- rep(NULL, n.rep)

        for (j in 1:n.rep){
          delta.data <- permutate.data[[j]][,valid.set]
          rep.bet[j] <- kmeans(delta.data, n.cluster, nstart = 5, iter.max = 50)$betweenss
        }

        gap[i] <- log(ori.between) - mean(log(rep.bet))
      }
      # the index of the position
      local.select <- grid.iter[which.max(gap)]
      # the index of the current iteration
      local.opt <- which.max(gap)
      gap.all[[iter]] <- gap
      if((final.pos - start.pos) / (maxnum - 1) <= 1){
        opt.n.noise <- local.select
        opt.gap <- max(gap)
      }
      else{
        if(local.select == start.pos){
          start.pos <- start.pos
          final.pos <- grid.iter[2]
        }
        else if(local.select == final.pos){
          start.pos <- grid.iter[length(grid.iter) - 1]
          final.pos <- final.pos
        }
        else{
          start.pos <- grid.iter[(local.opt - 1)]
          final.pos <- grid.iter[(local.opt + 1)]
        }
      }
    }
  }

  opt.results <- CardKMeans(dataset, n.cluster, opt.n.noise, kmeans_starts)
  opt.cluster.assign <- opt.results$cluster.assign
  opt.variables <- opt.results$variables

  results <- list(n.noisevar = opt.n.noise, signaling.set = opt.variables, cluster.assign = opt.cluster.assign, opt.cluster = n.cluster, org.data = dataset, gap = max(gap))
  class(results) <- "ckm"
  return(results)
}
