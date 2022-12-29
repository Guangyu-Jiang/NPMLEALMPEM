#' EM for solving the NPMLE with unknown supports
#' input: observations X \in R^{n*d}
#' SIGMA: n numbers of diagonal covararince matrices \in R^{d*n}
#' m: number of grid points
#' option$grid_initial: initial supports
#' @export
EM <- function(X, SIGMA, m, options){
  cat(sprintf('\n'))
  cat('----------------- EM algorithm--------------------')
  maxiter <- 100
  stoptol <- 1e-4
  n <- dim(X)[1]
  d <- dim(X)[2]
  ## uniform mixture proportions
  x <- 1/m*matrix(1, m, 1)##m*1
  if (m<n){
    supps <- X[sample(n,m),]##m*d
  }else if(m==n){
    supps <- X
  }
  if (!is.null(options$supps_initial)){
    supps <- options$supps_initial
  }
  if (!is.null(options$stoptol)){
    stoptol <- options$stoptol
  }
  Sigma <- t(SIGMA)
  inv_Sigma <- 1/Sigma##n*d
  inv_SigmaX <- inv_Sigma*X##n*d
  hist <- list()
  hist$obj <- rep(0, maxiter)
  for (k in 1:maxiter){
    L_list <- likelihood_matrix(X, supps, SIGMA, NULL, NULL)##n*m
    ##E step
    L <- L_list$L
    Lx <- L%*%x##n*1
    gamma_hat <- L*(matrix(1, n, 1)%*%t(x))/(Lx%*%matrix(1, 1, m))##n*m
    ##M step
    supps <- (t(gamma_hat)%*%inv_SigmaX)/(t(gamma_hat)%*%inv_Sigma)##m*d
    x <- matrix(colSums(gamma_hat), m, 1)/n##m*1
    obj <- sum(log(Lx))/n
    hist$obj[k] <- obj
    cat(sprintf('iter = %3d, log-likelihood = %5.8e \n', k, obj))
    if (k > 1){
      if(obj - obj_old  < stoptol){
      break
      }
    }
    obj_old <- obj

  }
  return(list(x = x, supps = supps, hist = hist, k = k))
}
