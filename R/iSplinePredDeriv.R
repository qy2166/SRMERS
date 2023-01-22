#' Calculate the derivatives of coefficients of a fixed effect regression spline (monotonicity).
#'
#' @export iSplinePredDeriv
#' @param t The knot sequence vector in exposure-outcome model.
#' @param theta The coefficient vector of I-spline bases in exposure-outcome model.
#' @param m The mediator value within min(t) and max(t).
#' @return The derivatives of beta.
iSplinePredDeriv <- function(t, theta, m){
  d <- rep(0, length(theta))
  for(k in 1:(length(theta)-1)){
    if(m < t[k+2] & m >= t[k+1]){
      fm = 0
      i = 1
      while(i < k){
        d[i] <- 1
        i = i + 1
      }
      d[k] <- (1 - (t[k+2]-m)^2/((t[k+2]-t[k+1])*(t[k+2]-t[k])))
      d[k+1] <- ((m-t[k+1])^2/((t[k+2]-t[k+1])*(t[k+3]-t[k+1])))
    }
  }
  names(d) <- names(theta)
  return(d)
}
