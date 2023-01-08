#' This function contains the main loop for Augmented Lagrange Multipler
#' @export
DualALM_main <- function(LL, parmain, x, y, u, v){
  tstart <- parmain$tstart
  stoptol <- parmain$stoptol
  stopop <- parmain$stopop
  printyes <- parmain$printyes
  maxiter <- parmain$maxiter
  approxL <- parmain$approxL
  approxRank <- parmain$approxRank
  sigma <- parmain$sigma
  m <- parmain$m
  n <- parmain$n
  stop <- 0
  sigmamax <- 1e7
  sigmamin <- 1e-8
  count_L <- 0
  count_LT <- 0
  ## initial objective values and feasibilities
  Lx <- LL$times(x)
  count_L <- count_L + 1
  obj <- c(sum(x) - sum(log(Lx))/n - 1, sum(log(u))/n)
  relgap <- abs(obj[1] - obj[2])/(1 + abs(obj[1]) + abs(obj[2]))
  Rp <- Lx - y
  normy <- norm(y, "2")
  primfeas <- max(norm(Rp, "2")/normy,norm(pmin(x,0), "2")/norm(x, "2"))
  LTv <- LL$trans(v)
  count_LT <- count_LT + 1
  Rd <- pmax(LTv - n, 0)
  normu <- norm(u, "2")
  dualfeas <- max(norm(Rd, "2")/n, norm(u - v, "2")/normu)
  maxfeas <- max(primfeas, dualfeas)
  eta <- norm(y - 1/v, "2")/normy
  if (printyes){
    cat(sprintf('\n (dimension: m = %d, n = %d, ', m, n))
    cat(sprintf('tol = %1.1e)\n', stoptol))
    cat('---------------------------------------------------')
    cat(sprintf('\n iter|  [pinfeas  dinfeas  complem]    relgap|'))
    cat(sprintf('    pobj          dobj         |  time sigma'))
    cat(sprintf('\n*********************************************'))
    cat(sprintf('*******************************************************'))
    cat(sprintf('\n%5.0f|  [%3.2e %3.2e %3.2e]   %3.2e| %- 8.6e %- 8.6e |',
                0, primfeas, dualfeas, eta, relgap, obj[1], obj[2]))
    cat(sprintf(' %5.1f  %3.2e', proc.time()[[3]]-tstart, sigma))
  }
  parNCG <- list()
  parNCG$tolconst <- 0.5
  parNCG$count_L <- count_L
  parNCG$count_LT <- count_LT
  parNCG$approxL <- approxL
  parNCG$approxRank <- approxRank
  parNCG$m <- m
  parNCG$n <- n
  maxitersub <- 20
  ssncgop <- list()
  ssncgop$tol <- stoptol
  ssncgop$printyes <- printyes

  ## main loop
  runhist<- list(primfeas = rep(NA, maxiter), dualfeas = rep(NA, maxiter),
         sigma = rep(NA, maxiter), primobj = rep(NA, maxiter),
         dualobj = rep(NA, maxiter), gap = rep(NA, maxiter),
         relgap = rep(NA, maxiter), ttime = rep(NA, maxiter),
         ttimessn = rep(NA, maxiter), itersub = rep(NA, maxiter),
         iterCG = rep(NA, maxiter), termination = NULL, iter = NULL
         )

  for (iter in 1:maxiter){
    parNCG$iter <- iter
    parNCG$sigma <- sigma
    if (dualfeas < 1e-5){
      maxitersub <- max(maxitersub, 35)
    }else if (dualfeas < 1e-3){
      maxitersub <- max(maxitersub, 30)
    }else if (dualfeas < 1e-1){
      maxitersub <- max(maxitersub, 30)
    }
    ssncgop$maxitersub <- maxitersub
    ## SSN
    tstart_ssn <- proc.time()[[3]]
    MLE_SSNCG_result <- MLE_SSNCG(LL, x, y, v, LTv, parNCG, ssncgop)
    x <- MLE_SSNCG_result$x
    y <- MLE_SSNCG_result$y
    u <- MLE_SSNCG_result$u
    v <- MLE_SSNCG_result$v
    Lx <- MLE_SSNCG_result$Lx
    LTv <- MLE_SSNCG_result$LTv
    parNCG <- MLE_SSNCG_result$par
    info_NCG <- MLE_SSNCG_result$info
    ttimessn <- proc.time()[[3]] - tstart_ssn
    if (info_NCG$breakyes < 0){
    parNCG$tolconst <- max(parNCG$tolconst/1.06, 1e-3)
    }
    ## compute KKT residual
    Rp <- Lx - y
    normy <- norm(y, "2")
    primfeas <- max(norm(Rp, "2")/normy, norm(min(x,0), "2")/norm(x, "2"))
    Rd <- pmax(LTv - n, 0)
    normu <- norm(u, "2")
    dualfeas <- max(norm(Rd, "2")/(n), norm(u - v, "2")/normu)
    maxfeas <- max(primfeas, dualfeas)
    eta <- norm(y - 1/v, "2")/normy
    ## compute objective values
    primobj <- sum(x) - sum(log(Lx))/n - 1
    dualobj <- sum(log(u))/n
    obj <- c(primobj, dualobj)
    gap <- primobj - dualobj
    relgap <- abs(gap)/(1 + abs(primobj) + abs(dualobj))
    if (stopop == 1 & maxfeas < stoptol & eta < stoptol){
      stop <- 1
    }else if (stopop == 2 & (eta < stoptol*10 | maxfeas < stoptol*10)){
      pkkt <- norm(x - pmax(x + LL$trans(1/(Lx))/n - 1, 0), "2")
      parNCG$count_LT <- parNCG$count_LT + 1
      if (pkkt < stoptol){
        stop <- 1
      }
    }else if (stopop == 3 & (eta < stoptol*10 | maxfeas < stoptol*10)){
      tmp <- LL$trans(1/(Lx))/n - 1
      parNCG$count_LT <- parNCG$count_LT + 1
      pkkt <- norm(x - pmax(x + tmp, 0), "2")
      pkkt2 <- max(tmp)
      if (max(pkkt2, pkkt) < stoptol){
        stop <- 1
      }
    }else if (stopop == 4){
      tmp <- LL$trans(1/(Lx))/n - 1
      parNCG$count_LT <- parNCG$count_LT + 1
      pkkt <- norm(x - max(x + tmp, 0), "2")
      if (pkkt < stoptol){
        stop <- 1
      }
    }
    ttime <- proc.time()[[3]] - tstart
    runhist$primfeas[iter] <- primfeas
    runhist$dualfeas[iter] <- dualfeas
    runhist$sigma[iter] <- sigma
    runhist$primobj[iter] <- primobj
    runhist$dualobj[iter] <- dualobj
    runhist$gap[iter] <- gap
    runhist$relgap[iter] <- relgap
    runhist$ttime[iter] <- ttime
    runhist$ttimessn[iter] <- ttimessn
    runhist$itersub[iter] <- info_NCG$itersub - 1
    runhist$iterCG[iter] <- info_NCG$tolCG
    if (printyes){
      cat(sprintf('\n%5.0f|  [%3.2e %3.2e %3.2e]   %3.2e| %- 8.6e %- 8.6e |',
                  iter,primfeas, dualfeas, eta, relgap, primobj, dualobj))
      cat(sprintf(' %5.1f  %3.2e', proc.time()[[3]]-tstart, sigma))
    }
    ## check termination
    if ((stop & iter > 5) | iter == maxiter){
      if (stop){
        termination <- 'converged'
        }else if (iter ==maxiter){
        termination <- 'maxiter reached'
        }
      runhist$termination <- termination
      runhist$iter <- iter
      obj[1] <- primobj
      obj[2] <- dualobj
      break
    }
    if (info_NCG$breakyes>=0) {
      sigma <- max(sigmamin, sigma/10)
      }else if (iter > 1){
      if(runhist$dualfeas[iter]/runhist$dualfeas[iter-1] > 0.6){
        if (sigma < 1e7 & primfeas < 100*stoptol){
        sigmascale <- 3
        }else {
        sigmascale = sqrt(3)
        }
        sigma <- min(sigmamax, sigma*sigmascale)
      }
    }
  }
  info <- list()
  info$maxfeas <- maxfeas
  info$eta <- eta
  info$iter <- iter
  info$relgap <- relgap
  info$ttime <- ttime
  info$termination <- termination
  info$msg <- termination
  info$count_L <- parNCG$count_L
  info$count_LT <- parNCG$count_LT
  return(list(obj = obj, x = x, y = y, u = u, v = v, info = info, runhist = runhist))
}
