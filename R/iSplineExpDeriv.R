#' Calculate the derivatives of coefficients of a fixed effect regression spline (monotonicity).
#'
#' @export iSplineExpDeriv
#' @param t The knot sequence vector in exposure-outcome model.
#' @param theta The coefficient vector of I-spline bases in exposure-outcome model.
#' @param sigma2 The residual standard deviation of exposure-mediator model.
#' @param gamma0 The coefficient of intercept of exposure-mediator model.
#' @param gamma1 The coefficient of exposure of exposure-mediator model.
#' @param gamma2 The coefficient vector of confounders of exposure-mediator model.
#' @param a The value of exposure (0/1).
#' @param c The values of confounders.
#' @return A list of derivatives of beta, gamma and sigma2^2.
iSplineExpDeriv <- function(t, theta, sigma2, gamma0, gamma1, gamma2, a, c){
  d1 <- vector()
  for(i in 1:length(theta)){
    f1 <- function(m){
      dtheta <- ((m-t[i])^2/((t[i+1]-t[i])*(t[i+2]-t[i])))
      dm <- (1/(sqrt(2*pi*sigma2^2)))*exp(-(m-(gamma0+gamma1*a+as.numeric(gamma2%*%c)))^2/(2*sigma2^2))
      return(dtheta*dm)
    }
    f2 <- function(m){
      dtheta <- (1 - (t[i+2]-m)^2/((t[i+2]-t[i+1])*(t[i+2]-t[i])))
      dm <- (1/(sqrt(2*pi*sigma2^2)))*exp(-(m-(gamma0+gamma1*a+as.numeric(gamma2%*%c)))^2/(2*sigma2^2))
      return(dtheta*dm)
    }
    f3 <- function(m){
      dm <- (1/(sqrt(2*pi*sigma2^2)))*exp(-(m-(gamma0+gamma1*a+as.numeric(gamma2%*%c)))^2/(2*sigma2^2))
      return(dm)
    }

    if(i == 1){
      d1[i] <- integrate(f2, t[i+1], t[i+2])$value
    }else if(i == length(theta)){
      d1[i] <- integrate(f1, t[i], t[i+1])$value
    }else{
      d1[i] <- integrate(f1, t[i], t[i+1])$value + integrate(f2, t[i+1], t[i+2])$value
    }

    j = i + 2
    while(j < (length(theta)+1)){
      d1[i] <- d1[i] + integrate(f3, t[j], t[j+1])$value
      j = j + 1
    }
  }
  names(d1) <- names(theta)

  d2 <- vector()
  dd2 <- c(-1, -a, -c)
  for(l in 1:length(dd2)){
    d2[l] <- 0
    for(k in 1:(length(theta)-1)){
      f <- function(m){
        gm = 0
        i = 1
        while(i < k){
          gm <- gm + theta[i]
          i = i + 1
        }
        gm <- gm + theta[k]*(1 - (t[k+2]-m)^2/((t[k+2]-t[k+1])*(t[k+2]-t[k]))) +
          theta[k+1]*((m-t[k+1])^2/((t[k+2]-t[k+1])*(t[k+3]-t[k+1])))
        dm <- (1/(sqrt(2*pi*sigma2^2)))*exp(-(m-(gamma0+gamma1*a+as.numeric(gamma2%*%c)))^2/(2*sigma2^2))
        ddm <- -(m-(gamma0+gamma1*a+as.numeric(gamma2%*%c)))/(sigma2^2)

        em <- gm*dm*ddm*dd2[l]

        return(em)
      }
      d2[l] <- d2[l] + integrate(f, t[k+1], t[k+2])$value
    }
  }
  names(d2) <- c(names(gamma0), names(gamma1), names(gamma2))

  d3 <- 0
  for(k in 1:(length(theta)-1)){
    f <- function(m){
      gm = 0
      i = 1
      while(i < k){
        gm <- gm + theta[i]
        i = i + 1
      }
      gm <- gm + theta[k]*(1 - (t[k+2]-m)^2/((t[k+2]-t[k+1])*(t[k+2]-t[k]))) +
        theta[k+1]*((m-t[k+1])^2/((t[k+2]-t[k+1])*(t[k+3]-t[k+1])))
      dm <- (1/(sqrt(2*pi*sigma2^2)))*exp(-(m-(gamma0+gamma1*a+as.numeric(gamma2%*%c)))^2/(2*sigma2^2))
      ddm <- -1/(2*sigma2^2) + (m-(gamma0+gamma1*a+as.numeric(gamma2%*%c)))^2/(2*(sigma2^2)^2)

      em <- gm*dm*ddm

      return(em)
    }
    d3 <- d3 + integrate(f, t[k+1], t[k+2])$value
  }
  names(d3) <- "sigma22"

  return(list(dbeta = d1,
              dgamma = d2,
              dsigma22 = d3))
}
