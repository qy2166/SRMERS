#' Calculate the predicted value of a fixed effect regression spline (monotonicity).
#'
#' @export iSplinePred
#' @param t The knot sequence vector in exposure-outcome model.
#' @param theta The coefficient vector of I-spline bases in exposure-outcome model.
#' @param m The mediator value within min(t) and max(t).
#' @return The predicted value of a fixed effect regression spline.
iSplinePred <- function(t, theta, m){
  if(m < min(t) | m >= max(t)){
    stop("mediator value is out of bounds")
  }else{
    for(k in 1:(length(theta)-1)){
      if(m < t[k+2] & m >= t[k+1]){
        fm = 0
        i = 1
        while(i < k){
          fm <- fm + theta[i]
          i = i + 1
        }
        fm <- fm + theta[k]*(1 - (t[k+2]-m)^2/((t[k+2]-t[k+1])*(t[k+2]-t[k]))) +
          theta[k+1]*((m-t[k+1])^2/((t[k+2]-t[k+1])*(t[k+3]-t[k+1])))
      }
    }
  }
  names(fm) <- "iSplinePrediction"
  return(fm)
}
