#' Calculate the predicted value of a fixed effect regression spline (convexity).
#'
#' @export cSplinePred
#' @param t The knot sequence vector in exposure-outcome model.
#' @param theta The coefficient vector of C-spline bases in exposure-outcome model.
#' @param m The mediator value within min(t) and max(t).
#' @return The predicted value of a fixed effect regression spline.
cSplinePred <- function(t, theta, m){
  if(m < min(t) | m >= max(t)){
    stop("mediator value is out of bounds")
  }else{
    thetaUpdate <- theta[-c(1)]
    for(k in 1:(length(thetaUpdate)-1)){
      if(m < t[k+2] & m >= t[k+1]){
        fm = 0
        i = 1
        while(i < k){
          fm <- fm + thetaUpdate[i]*(m - (t[i]+t[i+1]+t[i+2])/3)
          i = i + 1
        }
        fm <- fm + thetaUpdate[k]*(m - (t[k]+t[k+1]+t[k+2])/3 + (t[k+2]-m)^3/(3*(t[k+2]-t[k+1])*(t[k+2]-t[k]))) +
          thetaUpdate[k+1]*((m-t[k+1])^3/(3*(t[k+2]-t[k+1])*(t[k+3]-t[k+1])))
      }
    }
    fm <- fm + theta[c(1)]*m
  }
  names(fm) <- "cSplinePrediction"
  return(fm)
}
