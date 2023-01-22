#' Asymptotic null distribution of the test statistic
#'
#' @import nloptr
#' @import stats
#' @import MASS
#' @export distributionAsymp
#' @param nIter The number of iterations.
#' @param nBasis The number of bases.
#' @param sigmai The covariance matrix of coefficient estimates of I-Spline bases (monotonicity)
#' @param sigmac The covariance matrix of coefficient estimates of C-Spline bases (convexity)
#' @param testStatIncr The test statistic generated from iSplineFER/iSpineMER
#' @param testStatDecr The test statistic generated from iSplineFER/iSpineMER
#' @param testStatConv The test statistic generated from cSplineFER/cSpineMER
#' @param testStatConc The test statistic generated from cSplineFER/cSpineMER
#' @return A list of weights of beta distribution and p-values.
distributionAsymp <- function(nIter, nBasis, sigmai = NULL, sigmac = NULL, testStatIncr = NULL, testStatDecr = NULL, testStatConv = NULL, testStatConc = NULL){
  evalFun <- function(theta, thetaHat, sigma){
    objective <- t(theta - thetaHat) %*% solve(sigma) %*% (theta - thetaHat)
    return(objective)
  }
  x0 <- rep(0.1, nBasis)
  opts <- list("algorithm" = "NLOPT_LN_BOBYQA",
               "xtol_rel" = 1e-16)

  betaibar <- function(w, q, nBasis){
    fx <- 0
    for(i in 1:(nBasis+1)){
      fx <- fx + w[i]*pbeta(q, (i-1)/2, (nBasis-i+1)/2, lower.tail = F)
    }
    return(fx)
  }
  betacbar <- function(w, q, nBasis){
    fx <- 0
    for(i in 1:(nBasis+1)){
      fx <- fx + w[i]*pbeta(q, (i-1)/2, (nBasis-i+1)/2, lower.tail = F)
    }
    return(fx)
  }

  wi <- rep(0, nBasis+1)
  wc <- rep(0, nBasis+1)
  p_value_incr <- NULL
  p_value_decr <- NULL
  p_value_conv <- NULL
  p_value_conc <- NULL

  if(!is.null(sigmai) & !is.null(sigmac)){
    for(j in 1:nIter){
      betai <- mvrnorm(1, rep(0, nBasis), sigmai)

      betai_tilde <- nloptr(x0 = x0,
                            eval_f = evalFun,
                            lb = c(rep(0, nBasis)),
                            ub = c(rep(Inf, nBasis)),
                            opts = opts,
                            thetaHat = betai,
                            sigma = sigmai)$solution

      wi[sum(betai_tilde > 1e-8)+1] <- wi[sum(betai_tilde > 1e-8)+1] + 1

      if(j%%1000 == 0){
        print(j)
      }
    }
    if(!is.null(testStatIncr)){
      p_value_incr <- betaibar(wi/nIter, testStatIncr, nBasis)
    }
    if(!is.null(testStatDecr)){
      p_value_decr <- betaibar(wi/nIter, testStatDecr, nBasis)
    }

    for(j in 1:nIter){
      betac <- mvrnorm(1, rep(0, nBasis), sigmac)

      betac_tilde <- nloptr(x0 = x0,
                            eval_f = evalFun,
                            lb = c(rep(0, nBasis)),
                            ub = c(rep(Inf, nBasis)),
                            opts = opts,
                            thetaHat = betac,
                            sigma = sigmac)$solution

      wc[sum(betac_tilde > 1e-8)+1] <- wc[sum(betac_tilde > 1e-8)+1] + 1

      if(j%%1000 == 0){
        print(j)
      }
    }
    if(!is.null(testStatConv)){
      p_value_conv <- betaibar(wc/nIter, testStatConv, nBasis)
    }
    if(!is.null(testStatConc)){
      p_value_conc <- betaibar(wc/nIter, testStatConc, nBasis)
    }
  }
  else if(!is.null(sigmai) & is.null(sigmac)){
    for(j in 1:nIter){
      betai <- mvrnorm(1, rep(0, nBasis), sigmai)

      betai_tilde <- nloptr(x0 = x0,
                            eval_f = evalFun,
                            lb = c(rep(0, nBasis)),
                            ub = c(rep(Inf, nBasis)),
                            opts = opts,
                            thetaHat = betai,
                            sigma = sigmai)$solution

      wi[sum(betai_tilde > 1e-8)+1] <- wi[sum(betai_tilde > 1e-8)+1] + 1

      if(j%%1000 == 0){
        print(j)
      }
    }
    if(!is.null(testStatIncr)){
      p_value_incr <- betaibar(wi/nIter, testStatIncr, nBasis)
    }
    if(!is.null(testStatDecr)){
      p_value_decr <- betaibar(wi/nIter, testStatDecr, nBasis)
    }
  }
  else if(is.null(sigmai) & !is.null(sigmac)){
    for(j in 1:nIter){
      betac <- mvrnorm(1, rep(0, nBasis), sigmac)

      betac_tilde <- nloptr(x0 = x0,
                            eval_f = evalFun,
                            lb = c(rep(0, nBasis)),
                            ub = c(rep(Inf, nBasis)),
                            opts = opts,
                            thetaHat = betac,
                            sigma = sigmac)$solution

      wc[sum(betac_tilde > 1e-8)+1] <- wc[sum(betac_tilde > 1e-8)+1] + 1

      if(j%%1000 == 0){
        print(j)
      }
    }
    if(!is.null(testStatConv)){
      p_value_conv <- betaibar(wc/nIter, testStatConv, nBasis)
    }
    if(!is.null(testStatConc)){
      p_value_conc <- betaibar(wc/nIter, testStatConc, nBasis)
    }
  }
  else{
    stop("Please enter at least one sigma matrix")
  }

  return(list(iSplineWeight = wi,
              cSplineWeight = wc,
              PValueIncr = p_value_incr,
              PValueDecr = p_value_decr,
              PValueConv = p_value_conv,
              PValueConc = p_value_conc))
}
