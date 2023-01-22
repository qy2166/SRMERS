#' Calculate the expectation of a fixed effect regression spline (convexity).
#'
#' @export cSplineExp
#' @param t The knot sequence vector in exposure-outcome model.
#' @param theta The coefficient vector of C-spline bases in exposure-outcome model.
#' @param sigma2 The residual standard deviation of exposure-mediator model.
#' @param gamma0 The coefficient of intercept of exposure-mediator model.
#' @param gamma1 The coefficient of exposure of exposure-mediator model.
#' @param gamma2 The coefficient vector of confounders of exposure-mediator model.
#' @param a The value of exposure (0/1).
#' @param c The values of confounders.
#' @return The expectation of a fixed effect regression spline.
cSplineExp <- function(t, theta, sigma2, gamma0, gamma1, gamma2, a, c){
  thetaUpdate <- theta[-c(1)]
  expRes <- 0
  for(k in 1:(length(thetaUpdate)-1)){
    f <- function(m){
      gm = 0
      i = 1
      while(i < k){
        gm <- gm + thetaUpdate[i]*(m - (t[i]+t[i+1]+t[i+2])/3)
        i = i + 1
      }
      gm <- gm + thetaUpdate[k]*(m - (t[k]+t[k+1]+t[k+2])/3 + (t[k+2]-m)^3/(3*(t[k+2]-t[k+1])*(t[k+2]-t[k]))) +
        thetaUpdate[k+1]*((m-t[k+1])^3/(3*(t[k+2]-t[k+1])*(t[k+3]-t[k+1])))
      dm <- (1/(sqrt(2*pi*sigma2^2)))*exp(-(m-(gamma0+gamma1*a+as.numeric(gamma2%*%c)))^2/(2*sigma2^2))
      em <- gm*dm

      return(em)
    }
    expRes <- expRes + integrate(f, t[k+1], t[k+2])$value
  }
  expRes <- expRes + theta[c(1)]*(gamma0+gamma1*a+as.numeric(gamma2%*%c))
  names(expRes) <- "cSplineExpectation"
  return(expRes)
}
