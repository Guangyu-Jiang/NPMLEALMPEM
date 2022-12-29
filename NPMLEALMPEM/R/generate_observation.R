#' generate observations in R^2
#' Input: n # observations
#'        fig_option = 1 two concentric circles
#'                   = 2 triangle
#'                   = 3 digit 8
#'                   = 4 letter A
#'                   = 5 circle of radius 6
#'                   = 6 theta_i = 0
#'                   = 7 theta_i = 0, 6*e1, or 6*e2
#'                   = 8 theta_i = N(0,SIGMA);
#'                   = 9: theta_i has 6 support points [0  0 0 ... 0] [1  0 0 ... 0] [-1 0 0 ... 0]
#'                    [0  1 0 ... 0] [1  1 0 ... 0] [-1 1 0 ... 0]
#' sigma_option = 1, I2,
#'              = 2, SIGMA_i = diag(ai,bi),ai,bi ~ Uniform[1,3]
#' Output: observations X \in R^{n*d} (n observations in each row)
#'         theta \in R^{n*d} true signal
#'         SIGMA \in R^{d*n} if sigma_option = 2, SIGMA(:,:,i) = [ai,bi]
#' @export
generate_observation <- function(n, fig_option, sigma_option, d){
  if (is.null(d) | d < 2){
    d <- 2
  }
  if (is.null(sigma_option) | sigma_option ==1){
    SIGMA <- diag(d)
  }
  if(fig_option==1){
    r1 <- 2
    r2 <- 6
    n1 <- round(n/2)
    n2 <- n - n1
    theta <- matrix(0, n, d)
    t <- (2*pi)*runif(n1)
    theta[1:n1,1:2] <- r1*matrix(c(cos(t), sin(t)), n1, 2)
    t <- (2*pi)*runif(n2)
    theta[(n1 + 1):n,1:2] <- r2*matrix(c(cos(t), sin(t)), n2, 2)
    X <- t(apply(theta, 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
  }
  if(fig_option==2){
    n1 <- round(n/3)
    n2 <- n1
    n3 <- n - n1 - n2
    p1 <- c(-3,0)
    p2 <- c(0,6)
    p3 <- c(3,0)
    theta <- matrix(0, n, d)
    theta[,1:2] <- rbind(t(p1 + t(matrix(runif(n1),n1,1)%*%(p2-p1))),
                   t(p2 + t(matrix(runif(n2),n2,1)%*%(p3-p2))),
                   t(p3 + t(matrix(runif(n3),n3,1)%*%(p1-p3))))
    X <- t(apply(theta, 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
  }
  if(fig_option==3){
    r1 <- 3
    c1 <- c(0, 0)
    r2 <- 3
    c2 <- c(0, 6)
    n1 <- round(n/2)
    n2 <- n - n1
    theta <- matrix(0, n, d)
    t <- (2*pi)*runif(n1)
    theta[1:n1,1:2] <- c1 + r1*matrix(c(cos(t), sin(t)), n1, 2)
    t <- (2*pi)%*%runif(n2)
    theta[(n1 + 1):n,1:2] <- c2 + r2*matrix(c(cos(t), sin(t)), n2, 2)
    X <- t(apply(theta, 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
  }
  if(fig_option==4){
    n1 <- round(n/5)
    n2 <- n1
    n3 <- n1
    n4 <- n1
    n5 <- n - (n1 + n2 + n3 + n4)
    p1 <- c(-4, -6)
    p2 <- c(-2, 0)
    p3 <- c(0, 6)
    p4 <- c(2, 0)
    p5 <- c(4,-6)
    theta <- matrix(0, n, d)
    theta[,1:2] <- rbind(t(p1 + t(matrix(runif(n1), n1, 1)%*%(p2-p1))),
                   t(p2 + t(matrix(runif(n2), n2, 1)%*%(p3-p2))),
                   t(p3 + t(matrix(runif(n3), n3, 1)%*%(p4-p3))),
                   t(p4 + t(matrix(runif(n4), n4, 1)%*%(p5-p4))),
                   t(p4 + t(matrix(runif(n5), n5, 1)%*%(p2-p4))))
    X <- t(apply(theta, 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
  }
  if(fig_option==5){
    r <- 6
    t <- (2*pi)*runif(n)
    theta <- matrix(0, n, d)
    theta[,1:2] <- matrix(r*c(cos(t), sin(t)), n, 2)
    X <- t(apply(theta, 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
  }
  if(fig_option==6){
    theta <- matrix(0, n, d)
    X <- t(apply(theta, 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
  }
  if(fig_option==7){
    r <- 6
    atoms <- matrix(0, 3, d)
    atoms[2,1] <- r
    atoms[3,2] <- r
    theta <- matrix(0, n, d)
    xstar <- rep(1,3)/3## equal weights
    xstar2 <- cumsum(xstar)
    kk <- 0
    for (i in 1:3){
      kend <- round(n*xstar2[i])
      theta[(kk + 1):kend,1:2] <- atoms[rep(i,kend - kk),]
      kk <- kend
    }
    X <- t(apply(theta, 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
  }
  if(fig_option==8){
    theta <- matrix(0, n, d)
    theta[,1:2] = t(apply(matrix(0, n, d), 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
    X <- t(apply(theta, 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
  }
  if(fig_option==9){
    r = 6
    atoms <- matrix(0, 6, d)
    atoms[2,1] = r
    atoms[3,1] = -r
    atoms[4,2] = r
    atoms[5,1:2] = c(r, r)
    atoms[6,1:2] = c(-r, r)
    theta <- matrix(0, n, d)
    xstar <- rep(1, 6)/6## equal weights
    xstar2 <- cumsum(xstar)
    kk <- 0
    for (i in 1:6){
      kend <- round(n*xstar2[i])
      theta[(kk + 1):kend,1:2] <- atoms[rep(i,kend - kk),]
      kk <- kend
    }
    X <- t(apply(theta, 1, function(x){mvtnorm::rmvnorm(1, mean = x, sigma = SIGMA)}))
  }
  return(list(X = X, theta = theta, SIGMA = SIGMA))
}
