#' The function generates datasets that follow the typical K-means model with the option of including masking variables - variables that do not contribute to the clusters.
#' As a simplistic version, the current function restricts the means of all signaling variables to be equal within each cluster, while the variance to be equal across all variables and clusters.
#'
#' @param n.obs a positive integer indicating the number of observations
#' @param n.cluster a positive integer indicating the number of clusters
#' @param n.validvar an integer indicating the number of signaling variables
#' @param n.noisevar an integer indicating the number of masking variables
#' @param mu a vector of length \code{n.cluster} whole element indicates the mean value of each cluster. It could also be a number that indicates the distance of neighboring clusters
#' @param var a number indicates the variance of each variable
#' @param varsplit either 0 or 1 (default value is 0); when 1, the variance of half of the variables equal var/2
#' @return a list of two elements. The first is the generated dataset while the second is a vector of length \code{n.obs} contains the cluster assignment of each observations

#' @examples
#' ncluster <- 3
#' nobs <- 60
#' nnoisevar <- 100
#' nvalidvar <- 20
#' mu <- 1
#' var <- 1
#' sim.data <- DataGenCKM(nobs, ncluster, nvalidvar, nnoisevar, mu, var)
#' dataset <- sim.data[[1]]
#' cluster.assign <- sim.data[[2]]

DataGenCKM <- function(n.obs, n.cluster, n.validvar, n.noisevar, mu, var, varsplit = 0){

  if(length(mu) == 1){
    mu <- seq(from = -(n.cluster - 1) / 2 * mu, to = (n.cluster - 1) / 2 * mu, by = mu)
  }

  if (varsplit == 1){
    mu1 <- mu / 2
    n.validvar1 <- floor(n.validvar / 2)
  }

  n.var <- n.validvar + n.noisevar
  nk <- rep(NA, n.cluster)
  for(i in 1:(n.cluster-1)){
    nk[i] <- floor(n.obs/n.cluster)
  }
  nk[n.cluster] <- n.obs - floor(n.obs/n.cluster)*(n.cluster-1)

  ## create the score matrix
  partition.matrix <- matrix(0, nrow = n.obs, ncol = n.cluster)
  obs.set <- 1:n.obs

  for (i in 1:n.cluster){
    set.i <- sample(obs.set, nk[i])
    partition.matrix[set.i, i] <- 1
    obs.set <- obs.set[!(obs.set %in% set.i)]
  }

  cluster.assign <- sapply(1:n.obs, function(x) which(partition.matrix[x,] == 1))

  generate.matrix <- matrix(0, nrow = n.obs, ncol = n.var)

  if(varsplit == 0){
    for (j in 1:n.var){
      if (j < (n.validvar + 1)){
        for (i in 1:n.obs){
          i.cluster <- cluster.assign[i]
          generate.matrix[i,j] <- rnorm(1, mu[i.cluster], var)
        }
      }
      if(j == n.validvar){
        var.emp <- var(generate.matrix[,j])
      }
      if (j >= (n.validvar + 1)){
        generate.matrix[,j] <- rnorm(n.obs, 0, var)
      }
    }
  }

  if(varsplit == 1){
    for (j in 1:n.var){
      if(j < (n.validvar1 + 1)){
        for (i in 1:n.obs){
          i.cluster <- cluster.assign[i]
          generate.matrix[i,j] <- rnorm(1, mu[i.cluster], var)
        }
      }
      else if((j > n.validvar1) && j < (n.validvar + 1)){
        for (i in 1:n.obs){
          i.cluster <- cluster.assign[i]
          generate.matrix[i,j] <- rnorm(1, mu1[i.cluster], var)
        }
      }
      else{
        generate.matrix[,j] <- rnorm(n.obs, 0, var)
      }
    }
  }

  sd.matrix <- scale(generate.matrix)

  return(list(data = sd.matrix, cluster.assign = cluster.assign, mu = mu, ori.data = generate.matrix))
}
