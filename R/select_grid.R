#' select grid points
#' Input: observations X \in R^{n*d} (n observations in each row)
#'        grid_option = 1 grid points are chosen to be the data points
#'                = 2 grid points are taken as a random subsample of the data points
#'                = 3 use a uniform grid in a compact region (e.g.,rectangle in d=2) containing the data
#'                = 4 logspace of y, linspace of x
#' Output: grid points U \in R^{m*d} (m observations in each row)
#' @export
select_grid <- function(X, grid_option, m){
  U <- X
  n <- dim(X)[1]
  if (is.null(grid_option)){
    grid_option <- 1
  }
  if (grid_option == 1){
    U <- X
    m <- n
  }
  if (grid_option == 2){
    m <- min(m, n)
    U <- X[sample(n,m), ]
  }
  if (grid_option == 3){
    xmax <- max(X[, 1])
    xmin <- min(X[, 1])
    ymax <- max(X[, 2])
    ymin <- min(X[, 2])
    my <- round(sqrt(m))
    mx <- round(sqrt(m))
    xgrid <- seq(xmin, xmax, mx)
    ygrid = seq(ymin, ymax, my)
    meshgrid <- pracma::meshgrid(xgrid, ygrid)
    U <- cbind(c(meshgrid$X), c(meshgrid$Y))
  }
  if (grid_option == 4){
    xmax <- max(X[, 1])
    xmin <- min(X[, 1])
    ymax <- max(X[, 2])
    ymin <- min(X[, 2])
    my <- round(sqrt(m))
    mx <- round(sqrt(m))
    xgrid <- seq(xmin, xmax, mx)
    ygrid <- 10 ^ seq(from = log10(ymin), to = log10(ymax), length.out = my)
    #ygrid = pracma::logseq(ymin, ymax, (ymax-ymin)/my)
    meshgrid <- pracma::meshgrid(xgrid, ygrid)
    U <- cbind(c(meshgrid$X), c(meshgrid$Y))
  }
  return(list(U = U, m = m))
}

