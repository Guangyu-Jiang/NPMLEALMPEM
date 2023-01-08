#' generate the matrix L
#' Input: observations X \in R^{n*d}
#'        grid points U \in R^{m*d}
#'        SIGMA \in R^{1*d*n}
#'        normalizerows = 1, then max number of each row of output L is 1
#' Output: L \in R^{n*m}
#'         L_ij = mvtnorm(X[i,] - U[j,],SIGMA_i)
#'         rowmax[i] = max(L_ij,j=1,...,m)
#'         if normalizerows = 1, L[i,] = Lorig[i,:]/rowmax[i];
#'         elseif normalizerows = 0, L = Lorig.
#'         removeind denotes the index set of observations removed
#' @export
likelihood_matrix <- function(X, U, SIGMA, normalizerows, restrict_dist){
  if (is.null(normalizerows)){
    normalizerows <- 0
  }
  if (is.null(restrict_dist)){
    restrict_dist <- 0
  }
  n <- dim(X)[1]
  m <- dim(U)[1]
  L <- matrix(0, n, m)
  rowmax <- matrix(0, n, 1)
  if (restrict_dist){
    tiny <- 1e-9
  }else {
    tiny <- 1e-150
  }
  cnt <- 0
  removeind <- c()
  sz <- length(dim(SIGMA))
  for (i in 1:n){
    XI <- X[i, ]
    if (sz == 3){
      SIG <- SIGMA[1, , i]
    }else if (sz == 2){
      SIG <- SIGMA
    }else {
      SIG <- c()
    }
    tmp = mvtnorm::dmvnorm(x = t(XI-t(U)), sigma = diag(SIG[, i]))
    maxtmp <- max(tmp)
    if(maxtmp > tiny){
      cnt <- cnt + 1
      rowmax[cnt] <- maxtmp
      if(normalizerows){
        L[cnt,] <- pmax(tmp, tiny)/maxtmp
      }else{
        L[cnt,] <- pmax(tmp, tiny)
      }
    }else{
      removeind <- c(removeind, i)
    }
  }
  L <- L[1:cnt, ]
  rowmax=rowmax[1:cnt, ]
  return(list(L = L, rowmax = rowmax, removeind = removeind))
}
