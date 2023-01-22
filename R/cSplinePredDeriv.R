#' Calculate the derivatives of coefficients of a fixed effect regression spline (convexity).
#'
#' @export cSplinePredDeriv
#' @param t The knot sequence vector in exposure-outcome model.
#' @param theta The coefficient vector of C-spline bases in exposure-outcome model.
#' @param m The mediator value within min(t) and max(t).
#' @return The derivatives of beta.
cSplinePredDeriv <- function(t, theta, m){
  thetaUpdate <- theta[-c(1)]
  d <- rep(0, length(thetaUpdate))
  for(k in 1:(length(thetaUpdate)-1)){
    if(m < t[k+2] & m >= t[k+1]){
      fm = 0
      i = 1
      while(i < k){
        d[i] <- (m - (t[i]+t[i+1]+t[i+2])/3)
        i = i + 1
      }
      d[k] <- (m - (t[k]+t[k+1]+t[k+2])/3 + (t[k+2]-m)^3/(3*(t[k+2]-t[k+1])*(t[k+2]-t[k])))
      d[k+1] <- ((m-t[k+1])^3/(3*(t[k+2]-t[k+1])*(t[k+3]-t[k+1])))
    }
  }
  d <- c(d, m)
  names(d) <- c(names(theta)[-c(1)], names(theta)[c(1)])
  return(d)
}
