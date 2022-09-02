#' The function offers a summary of the results of CKM. The CKM analysis could either include the selection the number of clusters or without such an selection.
#' @param dataset a data frame that refers to the original data set
#' @param ckm.object an object of the S3 class "ckm", which could eitehr be the output of the function CKMSelVar (without the selection of number of clusters) or of the function CKMselAll (with the selection of nunber of clusters)
#' @return a list of five items that are corresponding to the selected or determine number of clusters. The first and second are integers that indicate the number of clusters and masking variables, respectively.
#' The third is a vector that contains the indices of signaling variables. The forth demonstrates the cluster assignment. The last one, a matrix, represent the centroids of the signaling variables (columns) for all clusters (rows).
#' @examples
#' ckm.sel.no <- CKMSelNo(dataset, ncluster, nnoisevar)
#' ckm.sel.var <- CKMSelVar(dataset, ncluster)
#' ckm.sel.all <- CKMSelAll(dataset, maxcluster)
#' summary(dataset, ckm.sel.no)
#' summary(dataset, ckm.sel.var)
#' summary(dataset, ckm.sel.all)
#'
#'
summary.ckm <- function(dataset, ckm.object){

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

  return(list(n.cluster = n.cluster, n.noisevar = n.noisevar, signaling.set = signaling.set, cluster.assign = cluster.assign, cluster.center = cluster.center))
}
