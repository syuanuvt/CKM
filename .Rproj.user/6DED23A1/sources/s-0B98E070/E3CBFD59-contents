#' The current function calculates the Cardinality KMeans (i.e. CKM) solution, without selecting the number of masking variables and clusters (they are assumed known)
#'
#' @param dataset the orginal dataset on which CKM operates
#' @param n.cluster the total number of clusters
#' @param n.noisevar the total number of irrelevant variables
#' @param num_starts_kmeans the number of starts for the conventional KM analysis (note that in CKM, the conventional KM operates in the lower dimensions). The default value is 20
#' @param recal whether a final step of KM is carried out to update cluster partitions (recommended when the number of signaling variable is small)
#' @return the function will return a ckm object that is the list of five elements. The first denotes the selected number of masking variables; the second includes all indicies of signaling variables; the third is a vector illustrating cluster assignment; the forth is the pre-determined or selected "optimal" number of clusters; the fifth is the original dataset.
#' @examples
#' ncluster <- 3
#' nnoisevar <- 100
#' ckm.sel.no <- CKMSelNo(dataset, ncluster, nnoisevar)

CKMSelNo <- function(dataset, n.cluster, n.noisevar, num_starts_kmeans = 20, recal = TRUE){
  ckm <- CardKMeans(dataset, n.cluster, n.noisevar, num_starts_kmeans, recal)
  opt.variables <- ckm$variables
  opt.cluster.assign <- ckm$cluster.assign
  results <- list(n.noisevar = n.noisevar, signaling.set = opt.variables, cluster.assign = opt.cluster.assign, opt.cluster = n.cluster)
  class(results) <- "ckm"
  return(results)
}
