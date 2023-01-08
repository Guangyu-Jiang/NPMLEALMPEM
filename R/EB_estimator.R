#' GMLEB estimator
#' Input the matrix L \in R^{n*m}
#'       vector x in R^m, x >= 0, sum xi = 1
#'       grid points U \in R^{m*d}
#' Output estimator theta_hat \in R^{n*d}
#' @export
EB_estimator <- function(L, x, U){
  Lx <- L %*% x##n*1
  n <- dim(L)[1]
  d <- dim(U)[2]
  theta_hat <- matrix(0, n, d)##n*d
  for (i in 1:n){
    theta_hat[i, ] <- (t(x)%*%(L[i, ]*U))/Lx[i]
  }
  return(theta_hat)
}
