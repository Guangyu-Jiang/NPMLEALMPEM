#' partial EM
#' input: observations X \in R^{n*d}
#' Sigma: n numbers of diagonal covaraince matrices \in R^{n*d}
#' m: number of grid points m<-n(PEM is mainly for d>=3)
#' mu: support set
#' option.grid_initial: initial supports
#' PEM for solving the NPMLE with unknown supports
#' input: observations X \in R^{n*d}
#' Sigma: n numbers of diagonal covaraince matrices \in R^{n*d}
#' m: number of grid points
#' option$grid_initial: initial supports
#' @export
PEM <- function(X, SIGMA, m, options){
  cat(' ----------------- partial EM algorithm--------------------')
  cat('\n')
  maxiter <- 100
  stoptol <- 1e-4
  printyes <- 1
  ##initialization of mu
  n <- dim(X)[1]##number of observations
  if (m <= n){
    supps <- X[sample(n, m), ]
    ##random subsample from observations X if m is no larger than n
  }
  if (!is.null(options$supps_initial)){
    supps <- options$supps_initial
  }
  if (!is.null(options$stoptol)){
    stoptol <- options$stoptol
  }
  if (!is.null(options$printyes)){
    printyes <- options$printyes
  }
  options$scaleL = 1
  options$approxL = 0
  options$printyes = printyes
  options$init_opt = 0
  options$stoptol = stoptol
  Sigma <- t(SIGMA)
  inv_Sigma <- 1/Sigma##n*d
  inv_SigmaX <- inv_Sigma*X##n*d
  ##PEM
  if (!is.null(options$rowmax)){
    options[['rowmax']] <- NULL
  }
  L <- likelihood_matrix(X, supps, SIGMA, NULL, NULL)$L##n*m
  hist <- list()
  hist$obj <- rep(0, maxiter)
  for (k in 1:maxiter){
    ##update x(weights)
    DualALM_result <- DualALM(L,options)
    obj <- DualALM_result$obj
    x <- DualALM_result$x##weights
    cat(sprintf("iter = %d, log-likelihood = %5.8e \n", k, -obj[1]))
    ##check convergence
    if (k>1){
      options$init_opt <- x
      if (obj[1]-obj_old<stoptol){
        break
      }
    }
    hist$obj[k] <- -obj[1]
    ##estimate supports
    posind <- x>0##m'
    sumposind <- sum(posind)
    xtmp <- x[posind]##positive weights m'
    Ltmp <- L[, posind]##L matrix corresponding to positive weights n*m'
    Lx <- Ltmp%*%xtmp
    gamma_hat <- Ltmp*(matrix(1, n, 1)%*%t(xtmp))/(Lx%*%matrix(1, 1, sumposind))##n*m
    supps_update <- (t(gamma_hat)%*%inv_SigmaX)/(t(gamma_hat)%*%inv_Sigma)##m*d
    supps[posind,] <- supps_update
    if (sumposind < m/3){
      L[,posind]=likelihood_matrix(X, supps_update, SIGMA, NULL, NULL)$L##update L
    }
    else {
      L=likelihood_matrix(X, supps, SIGMA, NULL, NULL)$L

    }
    obj_old <- obj[1]##update history
  }
  return (list(x = x, supps = supps, hist = hist, k = k))
}
