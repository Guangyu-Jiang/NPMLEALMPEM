#' Main function for Augmented Lagrange Multiplier
#' Input parameters:
#' L: a nonnegative n*m matrix
#' options: optional input parameter list containing several elements
#' options$stoptol is a tolerance parameter
#' options$printyes is a flag, =1 print details; =0 not print details
#' options$maxiter sets an upper bound on the number of iterations
#' options$alpha is the scaling parameter(by default, alpha = n)
#' Output parameters:
#' obj: primal and dual objective function values
#' x,y: primal variables
#' u,v: dual variables
#' info: information containing several fields
#' runhist: running history during the iterations, containing several fields
#' @export

#options<-list(stoptol = NULL,stopop = NULL,printyes = NULL,
#              maxiter = NULL,alpha = NULL,approxL = NULL,
#              approxRank = NULL, init_opt = NULL,rowmax = NULL)
#ROWMAX IS A VECTOR ABOUT L'S row max

# LL, parmain are created in this part
## LL includes LL$matrix;LL$U;LL$V  three matrixes and two functions LL$times() and LL$trans()


### function
DualALM <- function(L, options){
  stoptol <- 1e-6;   # error tolerance
  stopop <- 1;       # stopping option
  printyes <- 1;     # print details
  maxiter <- 100;    # maximum number of iterations
  sigma <- 100;      # penalty parameter
  scaleL <- 1;       # scale rows of L
  approxL <- 0;      # low rank approximation of L
  approxRank <- 30;
  init_opt <- 0;
  if (exists("options")) {
  # code here
    if(!is.null(options$stoptol)){
      stoptol <- options$stoptol
      }
    if(!is.null(options$stopop)){
      stopop <- options$stopop
      }
    if(!is.null(options$printyes)){
      printyes <- options$printyes
      }
    if(!is.null(options$maxiter)){
      maxiter <- options$maxiter
      }
    if(!is.null(options$sigma)){
      sigma <- options$sigma
      }
    if(!is.null(options$scaleL)){
      scaleL <- options$scaleL
      }
    if(!is.null(options$approxL)){
      approxL <- options$approxL
      }
    if(!is.null(options$approxL)){
      approxL <- options$approxL
      }
    if(!is.null(options$approxRank)){
      approxRank <- options$sapproxRank
      }
    if(!is.null(options$init_opt)){
      init_opt <- options$init_opt
    }
  }
  if (printyes) {
    cat(sprintf("\n*************************************************************************************"))
    cat(sprintf("\n ALM for the dual problem"))
    cat(sprintf("\n*************************************************************************************"))
    }
  ## tstart = clock; set a clock
  tstart <- proc.time()[[3]]
  n<-dim(L)[1]
  m<-dim(L)[2]

##function 1   scale L such that the maximal value of each row is 1
  if (scaleL) {
    if (exists("options") & (!is.null(options$rowmax))){
      s <- 1 / options$rowmax
    } else{
      s <- 1 / apply(L,1,max)
    }
    if (printyes){
      cat(sprintf('\n max/min scale = %3.1e/%3.1e', max(1./s), min(1./s)))
    }
    if (n > 1e6){
      for (i in 1:n){
        L[i, ] <- s[i] * L[i, ]
        }
      }else{
        L <- s * L
    }
  }else{
      s <- rep(1, n)
  }
  ##LL is created in this part
  LL <- list(matrix = L, U = NULL, V = NULL, times = NULL, trans = NULL)
  #construct LL as a list in this part

  ##function 2  low rank approximation of L
  if (approxL){
    approxRank <- ceiling(approxRank)
    approxSucceed <- 0
    t1 <- proc.time()[[3]]
    result <- svd(L, nu = approxRank, nv = approxRank)
    t2 <- proc.time()[[3]] - t1
    U <- result$u
    S <- diag(result$d)
    V <- result$v
    if (printyes) {
      cat(sprintf('\n approximate rank = %d',approxRank))
      cat(sprintf('\n partial svd(L) = %3.1f seconds',t2))
      cat(sprintf('\n eigenvalues <= %2.2e were truncated ', S[approxRank, approxRank]))
      print("\n ----------")
    }
    if (S[approxRank, approxRank] <= min(10 * stoptol, 1e-4)) {
      ii <- which(diag(S)< min(10 * stoptol, 1e-4))[1]
      if (!is.null(ii)) {
        U <- U[, 1:ii]
        S <- S[1:ii, 1:ii]
        V <- V[, 1:ii]
        approxRank <- ii
        if (printyes) {
          cat(sprintf('\n approximate rank = %d (further truncated)', approxRank))
          cat(sprintf('\n eigenvalues <= %2.2e were truncated ', S(approxRank, approxRank)))
          print("\n ----------")
        }
      }
      U <- U %*% S
      LL$U = U
      LL$V = V
      LL$times <- function(x) (U %*% (t(V) %*% x))
      LL$trans <- function(y) (V %*% (t(U) %*% y))
      approxSucceed <- 1
    }
    if (!approxSucceed & printyes) {
      cat(sprintf('\n numerical rank of L > %d, set approxL = 0', approxRank))
      print("\n ----------")
      approxL <- 0
      LL$times <- function(x) (L %*% x)
      LL$trans <- function(y) t((t(y) %*% L))
    }
  }else {
    LL$times <- function(x) (L %*% x)
    LL$trans <- function(y) t((t(y) %*% L))
  }
  ## function3  initialization
  if (init_opt == 0) {
    xnew <- matrix(1, m, 1)/m
    ynew <- rowSums(L) / m
    unew <- 1 / ynew
    vnew <- unew
  }else {
    xnew <- 0.5 * sigma * matrix(1, m, 1)
    ynew <- rowSums(L) / m
    unew <- 1 / ynew
    vnew <- matrix(2.2204e-16, n, 1)
  }
  Lx <- LL$times(xnew)
  tmp <- LL$trans(1 / Lx) / n - 1
  pkkt <- norm(xnew - pmax(xnew + tmp, 0), "2")
  if (pkkt < stoptol) {
    x <- xnew
    obj <- sum(x) + sum(log(s) - log(Lx)) / n - 1
    y <- numeric(0)
    u <- numeric(0)
    v <- numeric(0)
    info <- numeric(0)
    runhist <- numeric(0)
    cat(sprintf('\n Terminated at the initial point, primal KKT residual = %2.2e',pkkt))
    return(list(x = x, obj = obj, y = y, u = u, v = v, info = info, runhist = runhist))
  }
  ##stop at init
  ##Main algorithm, using the function ALM.MAIN
  parmain <- list(
    tstart = tstart,
    stoptol = stoptol,
    stopop = stopop,
    printyes = printyes,
    maxiter = maxiter,
    approxL = approxL,
    approxRank = approxRank,
    sigma = sigma,
    m = m,
    n = n
  )
  ## call Dual ALM_ main
  Main_return <- DualALM_main(LL, parmain, xnew, ynew, unew, vnew)
  xnew <- Main_return$x
  ynew <- Main_return$y
  unew <- Main_return$u
  vnew <- Main_return$v
  info_main <- Main_return$info
  runhist_main <- Main_return$runhist
  ttime <- proc.time()[[3]] - tstart
  iter <- info_main$iter
  msg <- info_main$msg

  if (iter == maxiter) {
    msg <- " maximum iteration reached"
  }
  x <- xnew
  y <- ynew / s
  u <- unew * s
  v <- vnew * s
  #compute original KKT residual
  Lx <- L %*% x
  Lxorg <- Lx / s
  Rp <- Lxorg - y
  normy <- norm(y, "2")
  primfeas = max(norm(Rp, "2") / norm(y, "2"), norm(pmin(x, 0), "2") / norm(x, "2"))
  Rd = pmax((t(vnew) %*% L) - n, 0)
  normu = norm(u, "2")
  dualfeas = max(norm(Rd, "2") / n, norm(u - v, "2") / normu)
  maxfeas = max(primfeas, dualfeas)
  eta = norm(y - 1 / v, "2") / norm(y, "2")
  #compute objective values
  primobj <- sum(x) + sum(log(s) - log(Lx))/n - 1
  dualobj <- sum(log(v))/n
  obj <- c(primobj, dualobj)
  gap <- primobj - dualobj
  relgap <- abs(gap)/(1 + abs(primobj) + abs(dualobj))
  tmp <- t((1/t(Lx))%*%L)/n -1
  pkkt <- norm(x - pmax(x + tmp, 0), "2")
  pkkt2 <- max(tmp)
  #Record infomation
  runhist <- runhist_main

  info <- list(relgap = relgap,
            iter = iter,
            itersub = sum(runhist$itersub),
            time = ttime,
            timessn = sum(runhist$ttimessn),
            eta = eta,
            obj = obj,
            maxfeas = maxfeas,
            kktres = max(maxfeas, eta),
            pkkt = pkkt,
            pkkt2 = pkkt2,
            sumlogLx = - sum(log(Lxorg)),
            count_L = info_main$count_L,
            count_LT = info_main$count_LT)
  if (printyes) {
    cat("\n****************************************\n")
    cat(sprintf('\n ALM          : %s', msg))
    cat(sprintf('\n iteration    : %d\n', iter))
    cat(sprintf(' L operator   : %d\n', info$count_L))
    cat(sprintf(' LT operator  : %d\n', info$count_LT))
    cat(sprintf(' time         : %3.2f\n', ttime))
    cat(sprintf(' prim_obj     : %4.8e\n', primobj))
    cat(sprintf(' dual_obj     : %4.8e\n', dualobj))
    cat(sprintf(' relgap       : %4.5e\n', relgap))
    cat(sprintf(' primfeas     : %3.2e\n', primfeas))
    cat(sprintf(' dualfeas     : %3.2e\n', dualfeas))
    cat(sprintf(' eta          : %3.2e\n', eta))
    cat(sprintf(' primalKKT    : %3.2e\n', pkkt))
    cat(sprintf(' primalKKT2   : %3.2e\n', pkkt2))
    cat(sprintf(' -sum(log(Lx)): %1.8e\n', info$sumlogLx))
    cat(sprintf(' sparsity     : %d\n', sum(x>0)))
  }
  return(list(x = x, obj = obj, y = y, u = u, v = v, info = info, runhist = runhist))
}






