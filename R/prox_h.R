#' proximal mapping of h(y) = (-1/n) \sum_j log(y_j), y \in R^n
#' prox_y = arg min_z { h(z) + (sigma/2) \|z - y\|^2 }
#' M_y = min_z { h(z) + (sigma/2) \|z - y\|^2 }
#' Input: y \in R^n, sigma \in R
#' Output:
#' prox_y \in R^n: the proximal point of y
#' M_y \in R: value of Moreau envelop at y
#' prox_prime \in R^n: diagonal vector of derivative of proximal maping
#' prox_prime_minus = 1 - prox_prime
#' @export

prox_h <- function(y, sigma){
  n <- length(y)
  tmp <- sqrt(y^2 + 4/(sigma*n))
  prox_y <- 0.5*(tmp + y)
  if (!is.na(which(prox_y<=0)[1])){
    cat(sprintf('\nlog is not defined in negative elements\n'))
    prox_y <- prox_y + 1e-30
  }
  M_y <- (sigma/2)%*%norm(prox_y - y, '2')^2 - (sum(log(prox_y)))/n
  tmp <- y/tmp
  prox_prime <- (1 + tmp)/2
  prox_prime_minus <- (1 - tmp)/2
  return(list(prox_y = prox_y, M_y = M_y, prox_prime = prox_prime,
              prox_prime_minus = prox_prime_minus))
}

