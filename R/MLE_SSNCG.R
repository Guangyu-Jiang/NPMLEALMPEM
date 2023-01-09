#' @export
MLE_SSNCG<- function(LL, x, y, v, LTv, par, options){
  eps <- 2.2204e-16
  printyes <- options$printyes
  maxitersub <- options$maxitersub
  tol <- options$tol
  breakyes <- 0
  maxitpsqmr <- 500
  precond <- 0
  stagnate_check_psqmr <- 0
  sigma <- par$sigma
  tiny <- 1e-10
  n <- par$n
  v1input <- v - (n/sigma)*y
  ## init run hist
  runhist<- list(priminf = rep(NA, maxitersub), dualinf = rep(NA, maxitersub),
                 phi = rep(NA, maxitersub),solve_ok = rep(NA, maxitersub),
                 psqmr = rep(NA, maxitersub),findstep = rep(NA, maxitersub))
  ## using function prox_h
  result_prox_h <- prox_h(v1input, sigma/(n^2))
  prox_v1 <- result_prox_h$prox_y
  M_v1 <- result_prox_h$M_y
  prox_v1_prime_m <- result_prox_h$prox_prime_minus
  v2input <- LTv/n + x/sigma - 1
  prox_v2 <- pmax(v2input, 0)
  Lprox_v2 <- LL$times(prox_v2)
  par$count_L <- par$count_L + 1
  phi <- -(M_v1 + (sigma/2)*norm(prox_v2, "2")^2)
  par$precond <- precond
  par$printyes <- printyes
  ##main Newton iteration
  for (itersub in 1:maxitersub) {
    tmp <- (sigma/n) * (v1input - prox_v1)
    Grad <- (tmp + sigma * Lprox_v2) / n
    normGrad <- norm(Grad, "2")
    priminf_sub <- normGrad / (norm(tmp, "2") / n)
    normu <- norm(prox_v1, "2")
    dualinf_sub <- max(norm(pmax(LTv - n, 0), "2") / (n),
                       norm(prox_v1 - v, "2") / normu)

    if(max(priminf_sub, dualinf_sub) < tol){
      tolsubconst <- 0.09
    }else{
      tolsubconst <- 0.005
    }

    tolsub <- max(min(1e-2, par$tolconst * dualinf_sub), tolsubconst * tol)
    runhist$priminf[itersub] <- priminf_sub
    runhist$dualinf[itersub] <- dualinf_sub
    runhist$phi[itersub] <- phi
    if (printyes) {
      cat(sprintf('\n      %2.0d  %- 11.10e  %3.2e %3.2e  %1.2e',
                  itersub, phi, priminf_sub, dualinf_sub, par$tolconst))
    }
    if (priminf_sub < tolsub & itersub > 1) {
      msg <- "good termination in subproblem:"
      if (printyes) {
        cat(sprintf('\n       %s  ',msg))
        cat(sprintf(' dualinf = %3.2e, normGrad = %3.2e, tolsub = %3.2e',
                    dualinf_sub, priminf_sub, tolsub))
      }
      u <- prox_v1
      x <- sigma * prox_v2
      Lx <- sigma * Lprox_v2
      y <- sigma/n * (prox_v1 - v1input)
      breakyes <- -1
      break
    }


    ##compute Newton direction
    ##precond = 0
    par$epsilon <- min(1e-3, 0.1 * normGrad)
    if (dualinf_sub > 1e-3 | itersub <= 5) {
      maxitpsqmr <- max(maxitpsqmr, 200)
    } else if (dualinf_sub > 1e-4) {
      maxitpsqmr <- max(maxitpsqmr, 300)
    } else if (dualinf_sub > 1e-5) {
      maxitpsqmr <- max(maxitpsqmr, 400)
    } else if (dualinf_sub > 5e-6) {
      maxitpsqmr <- max(maxitpsqmr, 500)
    }

    if (dualinf_sub > 1e-4) {
      stagnate_check_psqmr <- max(stagnate_check_psqmr, 20)
    } else {
      stagnate_check_psqmr <- max(stagnate_check_psqmr, 30)
    }
    if (itersub > 3 & dualinf_sub < 5e-5) {
      if(any(runhist$solve_ok[(itersub-3):(itersub -1)] <= -1)) {
        stagnate_check_psqmr <- max(stagnate_check_psqmr, 80)
      }
    }
    par$stagnate_check_psqmr <- stagnate_check_psqmr
    if (itersub > 1) {
      prim_ratio <- priminf_sub/runhist$priminf[itersub - 1]
      dual_ratio <- dualinf_sub/runhist$dualinf[itersub - 1]
    } else {
      prim_ratio <- 0
      dual_ratio <- 0
    }
    rhs <- -Grad
    if (par$iter < 2 & itersub < 5) {
      tolpsqmr <- min(1e-1, 0.01 * priminf_sub)
    } else {
      tolpsqmr <- min(1e-1, 0.001 * priminf_sub)
    }
    const2 <- 1
    if (itersub > 1 & (prim_ratio > 0.5 | priminf_sub > 0.1 * runhist$priminf[1])) {
      const2 <- 0.5 * const2
    }
    if (dual_ratio > 1.1) {
      const2 <- 0.5 * const2
    }
    tolpsqmr <- const2 * tolpsqmr
    par$tol <- tolpsqmr
    par$maxit <- maxitpsqmr
    par$minitpsqmr <- 5
    ##find Newton direction
    Linsolver_MLE_result <- Linsolver_MLE(rhs, LL, prox_v1_prime_m, v2input, par)
    dv <- Linsolver_MLE_result$dv
    resnrm <- Linsolver_MLE_result$resnrm
    solve_ok <- Linsolver_MLE_result$solve_ok
    par <- Linsolver_MLE_result$par
    iterpsqmr <- length(resnrm) - 1
    if (printyes) {
      cat(sprintf('| %3.1e %3.1e %3.0d %4d', par$tol, resnrm[length(resnrm)],
                  iterpsqmr, par$r))
      cat(sprintf(' %2.1f', const2))
    }
    ## line search
    if ((itersub <= 3) & (dualinf_sub > 1e-4) | (par$iter <= 3)) {
      stepop <- 1
    } else {
      stepop <- 2
    }
    steptol <- 10*1e-5
    LTdv <- LL$trans(dv)
    par$count_LT <- par$count_LT + 1
    Lprox_v2_new <- Lprox_v2
    findstep_return <- findstep(Grad, dv, LTdv, LL, phi, v1input, prox_v1,
                                prox_v1_prime_m, v2input, prox_v2,
                                Lprox_v2_new, steptol, stepop, par)
    phi <- findstep_return$phi
    v1input <- findstep_return$v1input
    prox_v1 <- findstep_return$prox_v1
    prox_v1_prime_m <- findstep_return$prox_v1_prime_m
    v2input <- findstep_return$v2input
    prox_v2 <- findstep_return$prox_v2
    Lprox_v2 <- findstep_return$Lprox_v2
    alp <- findstep_return$alp
    iterstep <- findstep_return$iter
    par <- findstep_return$par
    v <- v + alp*dv
    LTv <- LTv + alp*LTdv
    runhist$solve_ok[itersub] <- solve_ok
    runhist$psqmr[itersub] <- iterpsqmr
    runhist$findstep[itersub] <- iterstep
    if (alp < tiny) {
      breakyes <- 11
    }
    phi_ratio <- 1

    if (itersub > 1) {
      phi_ratio <- (phi - runhist$phi[itersub - 1])/(abs(phi) + eps)
    }
    if (printyes) {
      cat(sprintf(' %3.2e %2.0f', alp, iterstep))
      if (phi_ratio < 0) {
        #print('-')
      }
    }
    ##check for stagnation
    printsub <- printyes
    if (itersub > 4) {
      idx <- max(1, itersub - 3):itersub
      tmp <- runhist$priminf[idx]
      ratio <- min(tmp)/max(tmp)
      if ((any(runhist$solve_ok[idx] <= -1)) &
          (ratio > 0.9) &
          (min(runhist$psqmr[idx]) == max(runhist$psqmr[idx])) &
          (max(tmp) < 5*tol)) {
        #print('#')
        breakyes <- 1
      }
      const3 <- 0.7
      priminf_1half <- min(runhist$priminf[1:ceiling(itersub*const3)])
      priminf_2half <- min(runhist$priminf[ceiling(itersub*const3) + 1:itersub])
      priminf_best <- min(runhist$priminf[1:itersub - 1])
      priminf_ratio <- runhist$priminf[itersub]/runhist$priminf[itersub - 1]
      dualinf_ratio <- runhist$dualinf[itersub]/runhist$dualinf[itersub - 1]
      stagnate_idx <- which(runhist$solve_ok[1:itersub] <= -1)
      stagnate_count <- length(stagnate_idx)
      idx2 <- c(max(1, itersub - 7): itersub)
      if ((itersub >= 10) &
          (any(runhist$solve_ok[idx2] == -1)) &
          (priminf_best < 1e-2) &
          (dualinf_sub < 1e-3)) {
        tmp <- runhist$priminf[idx2]
        ratio <- min(tmp)/max(tmp)
        if (ratio > 0.5) {
          if (printsub) {
            #print('##')
          }
          breakyes <- 2
        }
      }
      if ((itersub >= 15) &
          (priminf_1half < min(2e-3, priminf_2half)) &
          (dualinf_sub < 0.8*runhist$dualinf[1]) &
          (dualinf_sub < 1e-3) &
          (stagnate_count >= 3)) {
        if (printsub) {
          #print('###')
        }
        breakyes <- 3
      }
      if ((itersub >= 15)&
         (priminf_ratio < 0.1) &
         (priminf_sub < 0.8*priminf_1half)&
         (dualinf_sub < min(1e-3, 2*priminf_sub))&
         ((priminf_sub < 2e-3)|((dualinf_sub < 1e-5) & (priminf_sub < 5e-3))) &
         (stagnate_count >= 3)){
        if (printsub) {
          #print(' $$')
        }
        breakyes <- 4
      }
      if (itersub >=10 & dualinf_sub > 5*min(runhist$dualinf[1:itersub]) &
          priminf_sub > 2*min(runhist$priminf[1:itersub])){
        if (printsub){
          #print('$$$')
        }
        breakyes <- 5
      }
      if (itersub >= 20) {
        dualinf_ratioall <- runhist$dualinf[2:itersub]/runhist$dualinf[1:itersub-1]
        idx <- which(dualinf_ratioall > 1)
        if (length(idx) >= 3) {
          dualinf_increment <- mean(dualinf_ratioall[idx])
          if (dualinf_increment > 1.25) {
            if (printsub) {
              #print('^^')
            }
            breakyes <- 6
          }
        }
      }
    }
  }
  if (itersub == maxitersub) {
    u <- prox_v1
    x <- sigma*prox_v2
    Lx <- sigma*Lprox_v2
    y <- sigma/n*(prox_v1 - v1input)
  }
  info <- list(tolCG = sum(runhist$psqmr[1:itersub]),
              breakyes= breakyes, itersub = itersub)
  return(list(x = x, y = y, u = u, v = v, Lx = Lx, LTv = LTv, par = par,
              runhist = runhist, info = info))
}
