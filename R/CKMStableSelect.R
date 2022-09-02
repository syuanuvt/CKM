CKMStableSelect <- function(dataset, n.cluster, n.noisevar, num_starts_kmeans = 20, recal = TRUE, nresample = 200, prop.resample = 1/2)
{
  num <- floor(n.noisevar * 1/2)
  var.seq <- floor(seq(from = num, to = n.noisevar, length.out = 10))
  test.grid <- expand.grid(nnoisevar = var.seq, number = 1:(nresample/10))
  sig.list <- matrix(nrow = nresample, ncol = (ncol(dataset)-num))
  for(i in 1:nresample){
    dataset.half <- dataset[sample(1:nrow(dataset), floor(nrow(dataset)*prop.resample)),]
    sig.list[i,1:(ncol(dataset)-test.grid[i,1])] <- CKMSelNo(dataset.half, n.cluster, test.grid[i,1], num_starts_kmeans = num_starts_kmeans)[[2]]
  }
  variables.rank <- table(sig.list)
  stable.set <- as.numeric(names(variables.rank)[sort.int(variables.rank,decreasing = TRUE,index.return = TRUE)$ix[1:(ncol(dataset)-n.noisevar)]])
  results <- list(stable.set = stable.set, rep.times = variables.rank[sort.int(variables.rank,decreasing = TRUE,index.return = TRUE)$ix])
  return(results)
}

