#' The function plots the results of CKM: the cluster centroids of all signaling variables
#' @param dataset a data frame that refers to the original data set
#' @param ckm.object an object of the S3 class "ckm", which could eitehr be the output of the function CKMSelVar (without the selection of number of clusters) or of the function CKMselAll (with the selection of nunber of clusters)
#' @param var.ind a vector of length J (i.e. the number of variables) that indicates the categories of variables
#' @return a plot with x axis representing various clusters and y axis referring to the centroids of signaling variables
#' @examples
#' ckm.sel.no <- CKMSelNo(dataset, ncluster, nnoisevar)
#' ckm.sel.var <- CKMSelVar(dataset, ncluster)
#' ckm.sel.all <- CKMSelAll(dataset, maxcluster)
#' plot(dataset, ckm.sel.no)
#' plot(dataset, ckm.sel.var)
#' plot(dataset, ckm.sel.all)
#'
#'
plot.ckm <- function(dataset, ckm.object, var.ind = NULL){

  n.noisevar <- ckm.object$n.noisevar
  signaling.set <- ckm.object$signaling.set
  cluster.assign <- ckm.object$cluster.assign
  n.cluster <- ckm.object$opt.cluster

  new.data <- dataset[, signaling.set]
  cluster.center <- matrix(nrow = n.cluster, ncol = length(signaling.set))
  colnames(cluster.center) <- as.character(signaling.set)

  for (i in 1:n.cluster){
    if(sum(cluster.assign == i) != 1){
      cluster.center[i,] <- apply(new.data[cluster.assign == i,], 2, mean)
    }
    if(sum(cluster.assign == i) == 1){
      cluster.center[i,] <- new.data[cluster.assign == i,]
    }
  }

  cluster.plot <- as.data.frame(t(cluster.center))
  for (i in 1:n.cluster){
    colnames(cluster.plot)[i] <- paste0("Cluster ", i)
  }

  if(is.null(var.ind)){

    print(GGally::ggparcoord(cluster.plot, columns = 1:n.cluster,  scale = "globalminmax",
             showPoints = TRUE, title = "Cluster Centroids", alphaLines = .3)+
      scale_color_viridis(discrete = TRUE)+
      ylab("Centroids") +
      xlab("Clusters"))
  }

  if(!is.null(var.ind)){
    cluster.plot$group <- var.ind[signaling.set]

    print(GGally::ggparcoord(cluster.plot, columns = 1:n.cluster, groupColumn = "group", scale = "globalminmax",
               showPoints = TRUE, title = "Cluster Centroids", alphaLines = .3)+
      scale_color_viridis(discrete = TRUE)+
      ylab("Centroids") +
      xlab("Clusters"))
  }
}
