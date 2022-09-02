#' The current function calculates the Cardinality KMeans (i.e. CKM) solution, without selecting the number of irrelevant variables
#' Note that the current function is an internal function
#'
#' @param dataset the orginal dataset on which CKM operates
#' @param n.cluster the total number of clusters
#' @param n.noisevar the total number of masking variables
#' @param num_starts_kmeans the number of starts for the conventional KM analysis. The default value is 10
#' @param recal whether a final step of KM is carried out to update cluster partitions (recommended when the number of signaling variable is small)
#' @return a list of three elements. The first is the vector that indicates the cluster assignment, the second is a vector that contains the
#' indices of all signaling variables, and the third is the total value of between-cluster sum of squares
#' @examples
#' ncluster <- 3
#' nnoisevar <- 100
#' CardKMeans(dataset, ncluster, nnoisevar)

CardKMeans <- function(dataset, n.cluster, n.noisevar, num_starts_kmeans = 10, recal = TRUE){
  global <- ComputeSCA(dataset, n.cluster, n.noisevar,  MAXITER = 100, stop_value = 1e-3, rational = 1)

  #############################################################################################
  selection.score <- global$score
  # the initial set of signaling variables
  variable.set <- (1:ncol(dataset))[!(apply(global$loading, 1, sum) == 0)]
  # a KM analysis conducted on the score matrix (i.e. lower dimensions)
  kmeans.results <- kmeans(selection.score, centers = n.cluster, nstart = num_starts_kmeans, iter.max = 100)$cluster

  ############################################################################################
  ## update the cluster assignment with the initial partition as a rational start
  new.data <- dataset[,!(apply(global$loading, 1, sum) == 0)]
  score.assign <- kmeans.results
  final.centers <- matrix(0, nrow = n.cluster, ncol = ncol(new.data))
  for (i in 1:n.cluster){
    if (sum(score.assign == i) != 1){
      final.centers[i, ] <- apply(new.data[score.assign == i, ], 2, mean)
    }
    else {
      final.centers[i, ] <- new.data[score.assign == i, ]
    }
  }
  kmeans.final <- kmeans(new.data, centers = n.cluster, iter.max = 100, nstart = num_starts_kmeans)
  kmeans.results <- kmeans.final$cluster

  ###########################################################################################
  ## further update the cluster assignment and the set of signaling variables in an iterative fashion
  conv <- 0
  converge <- 1e-6
  iter <- 0
  final.centers <- matrix(0, nrow = n.cluster, ncol = ncol(dataset))
  # calculate the initial cluster centers
  for (i in 1:n.cluster){
    if (sum(kmeans.results == i) != 1){
      final.centers[i, ] <- apply(dataset[kmeans.results == i, ], 2, mean)
    }
    else {
      final.centers[i, ] <- dataset[kmeans.results == i, ]
    }
  }
  minus.sum <- 0
  while(conv == 0){
    iter <- iter + 1
    distance <- rep(0, ncol(dataset))
    for (i in 1:nrow(dataset)){
      k <- kmeans.results[i]
      distance <- distance + dataset[i,] ^ 2 - (dataset[i,] - final.centers[k,]) ^ 2
    }
    variable.set <- sort(sort(distance, index.return = TRUE, decreasing = TRUE)$ix[1:(ncol(dataset) - n.noisevar)])
    obs.distance <- sum(distance[variable.set])

    for (i in 1:nrow(dataset)){
      distance.i <- sapply(1:n.cluster, function(x) (dataset[i,variable.set] - final.centers[x,variable.set])^2)
      if (is.matrix(distance.i)){
        kmeans.results[i] <- which(apply(distance.i, 2, mean) == min(apply(distance.i, 2, mean)))
      }
      else {
        kmeans.results[i] <- distance.i
      }
    }

    for (i in 1:n.cluster){
      if (sum(kmeans.results == i) != 1){
        final.centers[i, ] <- apply(dataset[kmeans.results == i, ], 2, mean)
      }
      else {
        final.centers[i, ] <- dataset[kmeans.results == i, ]
      }
    }

    if (iter == 1){
      new.loss <- obs.distance
      minus <- 1e9
    }
    else{
      old.loss <- new.loss
      new.loss <- obs.distance
      minus <- new.loss - old.loss
    }
    minus.sum <- c(minus.sum, minus)
    if (minus < converge)  conv <- 1
  }

  if(recal){
    ###########################################################################################
    ## calculate the final values that serve as the outputs
    kmeans.full <- kmeans(dataset[,variable.set], centers = n.cluster, nstart = num_starts_kmeans, iter.max = 50)
    if(kmeans.full$betweenss > new.loss){
      kmeans.results <- kmeans.full$cluster
      new.loss <- kmeans.full$betweenss
    }
  }

  return(list(cluster.assign = kmeans.results, variables = variable.set, betss = new.loss))
}
