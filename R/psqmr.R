#' psqmr:  preconditioned symmetric QMR with left (symmetric) preconditioner.
#' b = rhs vector.
#' resnrm = norm of qmr-generated residual vector b-Ax.
#' @export
psqmr <- function(matvecfname, b, par, x0, Ax0){
  N <- length(b)
  maxit <- max(5000, sqrt(length(b)))
  tol <- 1e-6 * norm(b, "2")
  stagnate_check <- 20
  miniter <- 0
  x0 <- matrix(0, N, 1)
  if (!is.null(par$maxit)){
    maxit <- par$maxit
  }
  if (!is.null(par$tol)){
    tol <- par$tol
  }
  if (!is.null(par$stagnate_check_psqmr)){
    stagnate_check <- par$stagnate_check_psqmr
  }
  if (!is.null(par$minitpsqmr)){
    miniter <- par$minitpsqmr
  }
  solve_ok <- 1
  printlevel <- 1
  x <- x0
  if (norm(x, "2") > 0){
    Ax0  <- matvecfname(x0)
    Aq <- Ax0
  }else{
    Aq <- matrix(0, N, 1)
  }
  r <- b - Aq
  err <- norm(r, '2')
  resnrm <- rep(0, maxit + 1)
  resnrm[1] <- err
  minres <- err
  q <- r
  tau_old  <- norm(q, '2')
  rho_old  <- t(r)%*%q
  theta_old <- 0
  d <- matrix(0, N, 1)
  res <- r
  Ad <- matrix(0, N, 1)
  ##main loop
  tiny <- -1e-30
  for (iter in 1:maxit){
    Aq <- matvecfname(q)
    sigma <- t(q)%*%Aq
    if (abs(sigma) < tiny){
      solve_ok <- 2
      if (printlevel){
        cat(sprintf('s1'))
      }
      break
    }else {
      alpha <- rho_old/sigma
      r <- r - c(alpha)*Aq
    }
    u <- r
    theta <- norm(u, "2")/tau_old
    c <- 1/sqrt(1 + theta^2)
    tau <- tau_old*theta*c
    gam <- (c^2*theta_old^2)
    eta <- (c^2*alpha)
    d <- gam*d + c(eta)*q
    x <- x + d
    ##----- stopping conditions ----
    Ad <- gam*Ad + c(eta)*Aq
    res <- res - Ad
    err <- norm(res, "2")
    resnrm[iter+1] <- err
    if (err < minres){
      minres <- err
    }
    if (err < tol & iter > miniter & t(b)%*%x > 0){
      break
    }
    if (iter > stagnate_check & iter > 10){
      ratio <- resnrm[iter-9:iter+1]/resnrm[iter-10:iter]
      if (min(ratio) > 0.997 & max(ratio) < 1.003){
        if (printlevel){
          cat(sprintf('s'))
        }
        solve_ok <- -1
        break
      }
    }
    if (abs(rho_old) < tiny){
      solve_ok <- 2
      cat(sprintf('s2'))
      break
      }else{
      rho <- t(r)%*%u
      beta <- rho/rho_old
      q <- u + c(beta)*q
      }
    rho_old <- rho
    tau_old <- tau
    theta_old <- theta
    }
  if (iter == maxit){
    solve_ok <- -2
    }
  Ax <- b - res
  return(list(x = x, Ax = Ax, resnrm = resnrm, solve_ok = solve_ok))
}


