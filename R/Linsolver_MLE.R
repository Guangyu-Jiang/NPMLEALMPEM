#' @export
Linsolver_MLE <- function(rhs, LL, prox_v1_prime_m, v2, par){
  n <- par$n
  m <- par$m
  sigma <- par$sigma
  eps <- 2.2204e-16
  J <- v2 > 0
  r <- sum(J)
  par$r <- r
  solveby <- 'pcg'  ## iterative solver cg
  if (n <= 5000){
    solveby <- 'pdirect'  ##direct solver
  }
  if (r < 2000 | (n > 5000 & r < 5000)){
    solveby = 'ddirect'  ## woodbury formula
  }
  if (r == 0){
    solveby <- 'none'
  }
  if (solveby == 'pdirect'){
    if (par$approxL){
      U <- LL$U
      V <- LL$V
      VJ <- V[J, ]
      LLT <- U%*%(t(VJ)%*%VJ)%*%t(U)
      for (i in 1:n){
        LLT[i, i] <- LLT[i, i] + prox_v1_prime_m[i]
      }
      cholLLT <- chol(LLT)
    }else{
      LJ <- LL$matrix[, J]
      LLT <- LJ%*%t(LJ)
      for (i in 1:n){
        LLT[i, i] <- LLT[i, i] + prox_v1_prime_m[i]
      }
      cholLLT <- chol(LLT)
    }
    dv <- mylinsysolve(cholLLT, rhs*(n^2/sigma))
    resnrm <- 0
    solve_ok <- 1
  }
  if (solveby == 'pcg'){
    if(par$approxL){
      U <- LL$U
      V <- LL$V
      VJ <- V[J, ]
      Afun <- function(v){
        return((prox_v1_prime_m*v + U%*%t(t(VJ%*%t(t(v)%*%U))%*%VJ))*(sigma/n^2))
      }
    }else {
      LJ <- LL$matrix[,J]
      Afun <- function(v){
        return((prox_v1_prime_m*v + LJ%*%(t(LJ)%*%v))*(sigma/n^2))
      }
    }
    result_psqmr <- psqmr(Afun, rhs, par, NULL, NULL)
    dv <- result_psqmr$x
    resnrm <- result_psqmr$resnrm
    solve_ok <- result_psqmr$solve_ok
  }
  if (solveby == 'ddirect'){
    par$approxL=0
    rhs <- rhs*(n^2/sigma)
    prox_v1_prime_m <- prox_v1_prime_m + eps
    if (par$approxL){
      U <- LL$U
      V <- LL$V
      VJ <- V[J,]
      rhstmp <- VJ%*%t(t(rhs/prox_v1_prime_m)%*%U)
      LTL <- VJ%*%(t(U)%*%(U/c(prox_v1_prime_m)))%*%t(VJ)
      for (i in 1:r){
        LTL[i,i] <- LTL[i,i] + 1
      }
      if (r<=1000){
        dv <- solve(LTL)%*%rhstmp
      }else {
        cholLTL <- chol(LTL)
        dv <- mylinsysolve(cholLTL, rhstmp)
      }
      dv <- (U%*%t(t(dv)%*%VJ))/prox_v1_prime_m
      dv <- rhs/prox_v1_prime_m - dv
    }else {
      if(r==m){
        LJ <- LL$matrix
      }
      else{
        LJ <-LL$matrix[ ,J]
      }
      LJ2 <- LJ/c(prox_v1_prime_m)
      rhstmp <- t(t(rhs)%*%LJ2)
      LTL <- t(LJ)%*%LJ2
      LTL <- diag(r) + LTL
      if (r <= 1000){
        dv <- solve(LTL)%*%rhstmp
      }else{
        cholLTL <- chol(LTL)
        dv <- mylinsysolve(cholLTL,rhstmp)
      }
      dv <- LJ2%*%dv
      dv <- rhs/prox_v1_prime_m - dv
    }
    resnrm <- 0
    solve_ok <- 1
  }
  if (solveby=='none'){
    dv <- rhs/prox_v1_prime_m*10000
    resnrm <- 0
    solve_ok <- 1
  }
  return(list(dv = dv, resnrm = resnrm, solve_ok = solve_ok, par = par))
}
