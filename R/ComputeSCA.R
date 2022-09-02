#' The function computes the solution of a special configuration of Sparce PCA (called SParse Variable PCA; SVPCA).
#' Instead of imposing zeros on individual elements of a loading matrix, all elements of a certain role (i.e. variables) are imposed to zero
#'
#' @param data the orginal dataset on which SVPCA operates
#' @param r number of components
#' @param n_sparvar number of variables that have all zero loadings
#' @param rational whether using rational starts (1, the default value) or random starts (0). The number of rational starts is fixed at 1. Otherwise, the number of starts is given by num_starts.
#' @param num_starts the number of random starts used when \code{rational == 0}. The default value is 10
#' @param MAXITER the maximum number of iterations. The default value is 1000
#' @param stop_value the convergence criteria. The default value if 1e-3
#' @return a list that contains three elements: the total loss value, the score matrices as well as the loading matrices
#' @examples
#' ncluster <- 3
#' nnoisevar <- 100
#' ComputeSCA(dataset, n.cluster, n.noisevar)

ComputeSCA<-function(data, r, n_sparvar, rational = 1, num_starts = 10, MAXITER = 1000, stop_value = 1e-3){

  n <- nrow(data)
  j <- ncol(data)

  zeros <- list()

  ## initialization of the loss function
  Loss_P <- 1e9
  min_Loss <- 1e9

  ## the initial values of P: SVD + a number randomly simulated from N(0, 1) distribution
  svd_data <- svds(data, r)
  initial_s <- svd_data$u
  #initial_p <- svd_data$v
  initial_p <- t(t(svd_data$v) * svd_data$d)
  global <- list()

  all_loss <- rep(0, num_starts)

  if(rational == 0){
    for (i in 1:num_starts){
      P <- initial_p + matrix(rnorm(j * r, 0, sqrt(sqrt(n))), nrow = j, ncol = r)#sqrt(sqrt(n))), nrow = j, ncol = r)
      #P <- initial_p + matrix(rnorm(j * r, 0, sqrt(n)), nrow = j, ncol = r)
      conv <- 0
      iter <- 0

      while (conv == 0)
      {
        #conditional estimation of T given P: closed form (see the book about LS estimates for detailed explanation)
        A <- t(data %*% P)
        svdA <- svd(A)
        T <- svdA$v %*% t(svdA$u)

        #update the current loss values
        DEV <- data - (T%*%t(P))
        DEVsq <- DEV^2
        fit <- sum(apply(DEVsq,2,sum))
        Loss_P <- fit

        ## update corresponding to each component (described in Gu & Van Deun, 2016)
        # conditional estimation of P given T
        P <- t(data) %*% T

        ## calculate the sum of squres of  each row, and impose zeros on the rows with smallest sum of squares
        row.ss <- apply(P^2, 1, sum)
        index.nsmallest <- sort(row.ss, index.return = TRUE)$ix[1:n_sparvar]
        P[index.nsmallest, ] <- 0

        #update the current loss values
        DEV <- data - (T%*%t(P))
        DEVsq <- DEV^2
        fit <- sum(apply(DEVsq,2,sum))
        update_loss <- (Loss_P - fit)
        Loss_P <- fit


        iter=iter+1

        ## stops
        if (iter > 1){
          if(update_loss < stop_value) conv <- 1
        }
        if(MAXITER==iter) conv <- 1
      }


      if (fit < min_Loss) {
        min_Loss <- fit
        global <- list(loss = min_Loss, score = T, loading = P)
      }

      all_loss[i] <- fit
      zeros[[i]] <- which(apply(P, 1, sum)==0)
    }
  }

  #one rational start
  if(rational == 1){
    P <- initial_p #+ matrix(rnorm(j * r, 0, sqrt(sqrt(n))), nrow = j, ncol = r)#sqrt(sqrt(n))), nrow = j, ncol = r)
    conv <- 0
    iter <- 0

    while (conv == 0)
    {
      #conditional estimation of T given P: closed form (see the book about LS estimates for detailed explanation)
      A <- t(data %*% P)
      svdA <- svd(A)
      T <- svdA$v %*% t(svdA$u)

      #update the current loss values
      DEV <- data - (T%*%t(P))
      DEVsq <- DEV^2
      fit <- sum(apply(DEVsq,2,sum))
      Loss_P <- fit

      ## update corresponding to each component (described in Gu & Van Deun, 2016)
      # conditional estimation of P given T
      P <- t(data) %*% T

      ## calculate the sum of squres of  each row, and impose zeros on the rows with smallest sum of squares
      row.ss <- apply(P^2, 1, sum)
      index.nsmallest <- sort(row.ss, index.return = TRUE)$ix[1:n_sparvar]
      P[index.nsmallest, ] <- 0

      #update the current loss values
      DEV <- data - (T%*%t(P))
      DEVsq <- DEV^2
      fit <- sum(apply(DEVsq,2,sum))
      update_loss <- (Loss_P - fit)
      Loss_P <- fit


      iter=iter+1

      ## stops
      if (iter > 1){
        if(update_loss < stop_value) conv <- 1
      }
      if(MAXITER==iter) conv <- 1
    }


    if (fit < min_Loss) {
      min_Loss <- fit
      global <- list(loss = min_Loss, score = T, loading = P)
    }
  }


  return(global)
}
