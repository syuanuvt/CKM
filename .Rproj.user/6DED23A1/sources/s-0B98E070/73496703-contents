#' The function selects the set of signaling variable, given the number of clusters
#'
#' @param dataset the orginal dataset on which CKM and its model selection procedure operates
#' @param minclust the minimal possible number of clusters
#' @param maxclust the maximal possible number of clusters
#' @param search the mode of selecting over the grid. "all" = selecting over each point of the grid; while it maximizes the accuracy, it is very slow with large number of variables.
#' "sub" = the "grid search with a zoom" strategy; while it is less accurate compared to searching the full grid, it is efficient even with large number of variables.
#' "dep" automatically adjust to one of the above two methods based on the number of variables. When # variables < 25, the search covers every possible value of the grid. This is also the default option.
#' @param maxnum the parameter is only useful when the "grid search with a zoom" strategy is applied. It restricts the maximal number of values searched over in any iteration. The default value is set at 10.
#' @param n.rep the number of permutated datasets when calculating the gap statistic
#' @param method different criterion exists as to determine the number of clusters based on the gap statistic; the users can try out these various options in the "method" argument. The default is "globalmax": selects the largest gap over all possible number of clusters (i.e. global maxima). Other options include "firstSEmax": select the "first" gap that falls within the range of the largest gap minus one SE (i.e. the one SE role); "firstmax": select the first largest gap (i.e., local maxima), "Tibs2001SEmax": the recommened guideline of Tibshirani, 2011 that takes the one-SE rule.
#' @param kmeans_starts the number of starts used in the kmeans algorithm. The default value is 10.
#' @param recal a logic that indicates whether a final step of KM is carried out to update cluster partitions (recommended when the number of signaling variable is small). The default value is TRUE.
#' @param sr a logic that indicates whether a scree ratio test has been added in the end (recommended when the number of clusters k is large (i.e., k > 15)). The default value is TRUE.
#' @param auto a logic that indicates whether the elbow point is determined automatically. Only useful when sr = TRUE. The default value is FALSE
#' @param ratio a numeric value that indicates the ratio used to identify the elbow point automatically. The elbow point is specified on which the difference in loss value of the current point is smaller than the ratio times the sum of the two previous differences in loss values
#' @return If auto = TRUE, the function will return a ckm object that is the list of four elements. The first denotes the selected number of masking variables; the second includes all indicies of signaling variables; the third is a vector illustrating cluster assignment; the forth is the pre-determined or selected "optimal" number of clusters.
#' If auto = FALSE, the function will plot loss values upon which the elbow point can be identified. With the found number of irrelevant variables, CKMSelNo can be used to determine the optimal solution.
#' @examples
#' mincluster <- 2
#' maxcluster <- 6
#' ckm.sel.all <- CKMSelAll(dataset, minclust, maxcluster)
#'
CKMSelVar <- function(dataset, n.cluster, search = "dep", maxnum = 10, n.rep = 20, kmeans_starts = 20, recal = TRUE, sr = TRUE, auto = FALSE, ratio = 1/3){

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
    for (i in search.grid) {
      method3.card <- try(CardKMeans(dataset, n.cluster, search.grid[i], kmeans_starts, recal))
      if(!inherits(method3.card,'try-error')){
        valid.set <- method3.card$variables
        ori.between <- method3.card$betss
        rep.bet <- rep(NULL, n.rep)
        for (j in 1:n.rep){
          delta.data <- permutate.data[[j]][,valid.set]
          rep.bet[j] <- kmeans(delta.data, n.cluster, nstart = 5, iter.max = 50)$betweenss
        }
        gap[i] <- log(ori.between) - mean(log(rep.bet))
      }
    }
    opt.n.noise <- which.max(gap)
  }

  if(search == "sub"){

    stop <- 0
    start.pos <- 1
    final.pos <- (n.var - 2 - floor(n.var / 500))
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
        method3.card <- try(CardKMeans(dataset, n.cluster, grid.iter[i], kmeans_starts, recal))
        if(!inherits(method3.card,'try-error')){
          valid.set <- method3.card$variables
          ori.between <- method3.card$betss
          rep.bet <- rep(NULL, n.rep)

          for (j in 1:n.rep){
            delta.data <- permutate.data[[j]][,valid.set]
            rep.bet[j] <- kmeans(delta.data, n.cluster, nstart = 5, iter.max = 50)$betweenss
          }

          gap[i] <- log(ori.between) - mean(log(rep.bet))
        }
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
    #opt.n.noise.old <- opt.n.noise
  }
    if(sr == TRUE){
      n.small <- opt.n.noise - 5
      n.big <- opt.n.noise + 5
      if(n.small < 1){
        n.small <- 1
      }
      if(n.big >= (ncol(dataset)-2)){
        n.big <- ncol(dataset)-2
      }
      alt.n.noise <- n.big:n.small
      alt.leng <- length(alt.n.noise)
      ori <- rep(NA,alt.leng)
      ind <- 0
      opt.set <- list()
      for (i in alt.n.noise) {
        ind <- ind + 1
        opt.set[[ind]] <- CardKMeans(dataset, n.cluster, i, n.rep, recal)
        ori[ind] <- opt.set[[ind]]$betss
      }
      ori.dif <- diff(ori)
      if(auto == TRUE){
        sign <- 0
        ind <- 2
        while(sign == 0){
          ind <- ind +1
          if(ori.dif[ind] < (ori.dif[ind-2]+ori.dif[ind-1])*ratio){
            sign <- 1
          }
          if(sign == 0 & ind == length(ori.dif)){
            sign <- 1
            ind <- which(alt.n.noise == opt.n.noise)
          }
        }
        index <- ind
        opt.n.noise <- alt.n.noise[index]
        opt.cluster.assign <- opt.set[[index]]$cluster.assign
        opt.variables <- opt.set[[index]]$variables
        results <- list(n.noisevar = opt.n.noise, signaling.set = opt.variables, cluster.assign = opt.cluster.assign, opt.cluster = n.cluster)
        class(results) <- "ckm"
        return(results)
      }
      if(auto == FALSE){
        plot.loss <- as.data.frame(cbind(rev(alt.n.noise), ori))
        names(plot.loss) <- c("nnoise","loss")
        print(ggplot(data = plot.loss, aes(x=nnoise, y=loss)) +
          geom_point()+
          scale_x_continuous(limits = c(min(alt.n.noise),max(alt.n.noise)), breaks = alt.n.noise) +
          theme(axis.title=element_text(size=12,face="bold"),
          panel.border= element_blank(),
        panel.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        legend.text = element_text(hjust = 1),
        legend.title = element_blank(),
        legend.key.width = unit(.3, "cm"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        axis.text.y = element_text( size = 10, color = "black"),
        axis.text.x = element_text( size = 10, color = "black", angle = 90, hjust = 1)))
      }
    }

    if(sr == FALSE){
      results <- CKMSelNo(dataset, n.cluster, opt.n.noise, kmeans_starts, recal)
      return(results)
    }
  }
