#' Calculate the CDE, NDE and NIE.
#'
#' @import splines2
#' @import coneproj
#' @import Matrix
#' @export SRSplineMed
#' @param data A data frame.
#' @param nBasis The number of bases.
#' @param exposure The name of the exposure (must be a binary variable).
#' @param mediator The name of the mediator (must be a continuous variable).
#' @param outcome The name of the outcome (must be a continuous variable).
#' @param confounderVec The name vector of the confounders.
#' @param shapeExp The shape of mediator in exposure group ("increasing", "decreasing", "convex", or "concave").
#' @param shapeNonExp The shape of mediator in non-exposure group ("increasing", "decreasing", "convex", or "concave").
#' @param mValue The controlled mediator value for CDE estimation.
#' @param varAsymp Whether to output the asymptotic variance (T/F)
#' @return A list of exposure-outcome model, exposure-mediator model, knot sequence, coefficient vector of exposure spline, coefficient vector of non-exposure spline, residuals, sds and coefficients, CDE, NDE and NIE and their asymptotic variances.
SRSplineMed <- function(data, nBasis, exposure, mediator, outcome, confounderVec, shapeExp, shapeNonExp, mValue, varAsymp = FALSE){
  m = data[, mediator]
  y = data[, outcome]
  a = data[, exposure]

  ## exposure > mediator > outcome
  knots <- quantile(m, seq(0, 1, 1/(nBasis - 1)), na.rm = T)
  knots <- knots[2:(length(knots)-1)]

  ## w
  w = data[, c(confounderVec, exposure)]

  ## z0 and z1
  if(shapeExp == "increasing"){
    z1 = (iSpline(m, knots = knots, degree = 1, intercept = T))*a
    z1m = (iSpline(m, knots = knots, degree = 1, intercept = T))*a
  }else if(shapeExp == "decreasing"){
    z1 = (1-iSpline(m, knots = knots, degree = 1, intercept = T))*a
    z1m = (iSpline(m, knots = knots, degree = 1, intercept = T))*a
  }else if(shapeExp == "convex"){
    z1 = (cSpline(m, knots = knots, degree = 1, intercept = T, scale = T))*a
    z1m = (cSpline(m, knots = knots, degree = 1, intercept = T, scale = F))*a
  }else if(shapeExp == "concave"){
    z1 = (1-cSpline(m, knots = knots, degree = 1, intercept = T, scale = T))*a
    z1m = (cSpline(m, knots = knots, degree = 1, intercept = T, scale = F))*a
  }else{
    stop("Please enter the correct value")
  }
  colnames(z1) <- paste0("splineExp", 1:ncol(z1))
  colnames(z1m) <- paste0("splineExp", 1:ncol(z1m))

  if(shapeNonExp == "increasing"){
    z0 = (iSpline(m, knots = knots, degree = 1, intercept = T))*(1-a)
    z0m = (iSpline(m, knots = knots, degree = 1, intercept = T))*(1-a)
  }else if(shapeNonExp == "decreasing"){
    z0 = (1-iSpline(m, knots = knots, degree = 1, intercept = T))*(1-a)
    z0m = (iSpline(m, knots = knots, degree = 1, intercept = T))*(1-a)
  }else if(shapeNonExp == "convex"){
    z0 = (cSpline(m, knots = knots, degree = 1, intercept = T, scale = T))*(1-a)
    z0m = (cSpline(m, knots = knots, degree = 1, intercept = T, scale = F))*(1-a)
  }else if(shapeNonExp == "concave"){
    z0 = (1-cSpline(m, knots = knots, degree = 1, intercept = T, scale = T))*(1-a)
    z0m = (cSpline(m, knots = knots, degree = 1, intercept = T, scale = F))*(1-a)
  }else{
    stop("Please enter the correct value")
  }
  colnames(z0) <- paste0("splineNonExp", 1:ncol(z0))
  colnames(z0m) <- paste0("splineNonExp", 1:ncol(z0m))

  z = cbind(z1, z0)
  zm = cbind(z1m, z0m)

  ## w0
  if(shapeNonExp %in% c("increasing", "decreasing") & shapeExp %in% c("increasing", "decreasing")){
    w0 = matrix(rep(1, nrow(data)), ncol = 1)
    colnames(w0) <- "Intercept"
  }else if(shapeNonExp %in% c("increasing", "decreasing") & shapeExp %in% c("convex", "concave")){
    w0 = cbind(1, m*a)
    colnames(w0) <- c("Intercept", paste0("splineExp0"))
  }else if(shapeNonExp %in% c("convex", "concave") & shapeExp %in% c("increasing", "decreasing")){
    w0 = cbind(1, m*(1-a))
    colnames(w0) <- c("Intercept", paste0("splineNonExp0"))
  }else if(shapeNonExp %in% c("convex", "concave") & shapeExp %in% c("convex", "concave")){
    w0 = cbind(1, m*a, m*(1-a))
    colnames(w0) <- c("Intercept", paste0("splineExp0"), paste0("splineNonExp0"))
  }else{
    stop("Please enter the correct value")
  }

  ## cone projection
  v = as.matrix(cbind(w0, w))

  # pv = v%*%solve(t(v)%*%v)%*%t(v)
  # delta = (diag(1, nrow = nrow(pv)) - pv) %*% z

  pvz_temp1 = v%*%solve(t(v)%*%v)
  pvz_temp2 = t(v)%*%z
  pvz = pvz_temp1 %*% pvz_temp2
  delta = z - pvz

  cone <- coneB(y, delta)

  ## linear regression
  zn <- zm[, which(cone$coefs != 0)]

  newdata1 <- cbind(y, zn, v[, -1])
  newdata1 <- data.frame(newdata1)
  colnames(newdata1)[1] <- "y"
  lm1 <- lm(y ~ ., data = newdata1)
  sigma1 <- sigma(lm1)
  residual1 <- residuals(lm1)

  coefExp <- lm1$coefficients[grep("splineExp", names(lm1$coefficients))]
  if(shapeExp %in% c("convex", "concave")){
    coefExpName <- paste0("splineExp", 0:nBasis)
  }else if(shapeExp %in% c("increasing", "decreasing")){
    coefExpName <- paste0("splineExp", 1:nBasis)
  }else{
    stop("Please enter the correct value")
  }
  coefExpUpdate <- c(coefExp, rep(0, length(setdiff(coefExpName, names(coefExp)))))
  names(coefExpUpdate) <- c(names(coefExp), setdiff(coefExpName, names(coefExp)))
  coefExpUpdate <- coefExpUpdate[order(factor(names(coefExpUpdate), levels = coefExpName))]

  coefNonExp <- lm1$coefficients[grep("splineNonExp", names(lm1$coefficients))]
  if(shapeNonExp %in% c("convex", "concave")){
    coefNonExpName <- paste0("splineNonExp", 0:nBasis)
  }else if(shapeNonExp %in% c("increasing", "decreasing")){
    coefNonExpName <- paste0("splineNonExp", 1:nBasis)
  }else{
    stop("Please enter the correct value")
  }
  coefNonExpUpdate <- c(coefNonExp, rep(0, length(setdiff(coefNonExpName, names(coefNonExp)))))
  names(coefNonExpUpdate) <- c(names(coefNonExp), setdiff(coefNonExpName, names(coefNonExp)))
  coefNonExpUpdate <- coefNonExpUpdate[order(factor(names(coefNonExpUpdate), levels = coefNonExpName))]

  ## exposure > mediator
  newdata2 <- data[, c(mediator, exposure, confounderVec)]
  colnames(newdata2)[1] <- "m"
  lm2 <- lm(m ~., data = newdata2)
  sigma2 <- sigma(lm2)
  residual2 <- residuals(lm2)
  gamma0 <- lm2$coefficient["(Intercept)"]
  gamma1 <- lm2$coefficient[exposure]
  gamma2 <- lm2$coefficient[confounderVec]

  t = c(min(m), quantile(m, seq(0, 1, 1/(nBasis - 1)), na.rm = T), max(m))

  c = apply(lm2$model, 2, mean)[confounderVec]

  if(shapeExp %in% c("convex", "concave")){
    expectExp0 <- cSplineExp(t, coefExpUpdate, sigma2, gamma0, gamma1, gamma2, 0, c)
    expectExp1 <- cSplineExp(t, coefExpUpdate, sigma2, gamma0, gamma1, gamma2, 1, c)
    predExp <- cSplinePred(t, coefExpUpdate, mValue)
  }else if(shapeExp %in% c("increasing", "decreasing")){
    expectExp0 <- iSplineExp(t, coefExpUpdate, sigma2, gamma0, gamma1, gamma2, 0, c)
    expectExp1 <- iSplineExp(t, coefExpUpdate, sigma2, gamma0, gamma1, gamma2, 1, c)
    predExp <- iSplinePred(t, coefExpUpdate, mValue)
  }else{
    stop("Please enter the correct value")
  }

  if(shapeNonExp %in% c("convex", "concave")){
    expectNonExp0 <- cSplineExp(t, coefNonExpUpdate, sigma2, gamma0, gamma1, gamma2, 0, c)
    expectNonExp1 <- cSplineExp(t, coefNonExpUpdate, sigma2, gamma0, gamma1, gamma2, 1, c)
    predNonExp <- cSplinePred(t, coefNonExpUpdate, mValue)
  }else if(shapeNonExp %in% c("increasing", "decreasing")){
    expectNonExp0 <- iSplineExp(t, coefNonExpUpdate, sigma2, gamma0, gamma1, gamma2, 0, c)
    expectNonExp1 <- iSplineExp(t, coefNonExpUpdate, sigma2, gamma0, gamma1, gamma2, 1, c)
    predNonExp <- iSplinePred(t, coefNonExpUpdate, mValue)
  }else{
    stop("Please enter the correct value")
  }
  CDE = lm1$coefficient[exposure] + predExp - predNonExp
  names(CDE) <- "CDE"
  NDE = lm1$coefficient[exposure] + expectExp0 - expectNonExp0
  names(NDE) <- "NDE"
  NIE = expectExp1 - expectExp0
  names(NIE) <- "NIE"

  if(varAsymp == F){
    CDEv <- NA
    names(CDEv) <- "CDE (var)"
    NDEv <- NA
    names(NDEv) <- "NDE (var)"
    NIEv <- NA
    names(NIEv) <- "NIE (var)"
  }else{
    if(shapeExp %in% c("convex", "concave")){
      derivExpectExp0 <- cSplineExpDeriv(t, coefExpUpdate, sigma2, gamma0, gamma1, gamma2, 0, c)
      derivExpectExp1 <- cSplineExpDeriv(t, coefExpUpdate, sigma2, gamma0, gamma1, gamma2, 1, c)
      derivPredExp <- cSplinePredDeriv(t, coefExpUpdate, mValue)
    }else if(shapeExp %in% c("increasing", "decreasing")){
      derivExpectExp0 <- iSplineExpDeriv(t, coefExpUpdate, sigma2, gamma0, gamma1, gamma2, 0, c)
      derivExpectExp1 <- iSplineExpDeriv(t, coefExpUpdate, sigma2, gamma0, gamma1, gamma2, 1, c)
      derivPredExp <- iSplinePredDeriv(t, coefExpUpdate, mValue)
    }else{
      stop("Please enter the correct value")
    }

    if(shapeNonExp %in% c("convex", "concave")){
      derivExpectNonExp0 <- cSplineExpDeriv(t, coefNonExpUpdate, sigma2, gamma0, gamma1, gamma2, 0, c)
      derivExpectNonExp1 <- cSplineExpDeriv(t, coefNonExpUpdate, sigma2, gamma0, gamma1, gamma2, 1, c)
      derivPredNonExp <- cSplinePredDeriv(t, coefNonExpUpdate, mValue)
    }else if(shapeNonExp %in% c("increasing", "decreasing")){
      derivExpectNonExp0 <- iSplineExpDeriv(t, coefNonExpUpdate, sigma2, gamma0, gamma1, gamma2, 0, c)
      derivExpectNonExp1 <- iSplineExpDeriv(t, coefNonExpUpdate, sigma2, gamma0, gamma1, gamma2, 1, c)
      derivPredNonExp <- iSplinePredDeriv(t, coefNonExpUpdate, mValue)
    }else{
      stop("Please enter the correct value")
    }
    derivPredExp <- derivPredExp[names(derivPredExp) %in% names(coefExp)]
    derivPredNonExp <- derivPredNonExp[names(derivPredNonExp) %in% names(coefNonExp)]
    derivExpectExp0$dbeta <- derivExpectExp0$dbeta[names(derivExpectExp0$dbeta) %in% names(coefExp)]
    derivExpectNonExp0$dbeta <- derivExpectNonExp0$dbeta[names(derivExpectNonExp0$dbeta) %in% names(coefNonExp)]
    derivExpectExp1$dbeta <- derivExpectExp1$dbeta[names(derivExpectExp1$dbeta) %in% names(coefExp)]
    derivExpectNonExp1$dbeta <- derivExpectNonExp1$dbeta[names(derivExpectNonExp1$dbeta) %in% names(coefNonExp)]

    lm1vcov <- vcov(lm1)
    lm2vcov <- vcov(lm2)

    v11 <- lm1vcov[c(grep(exposure, colnames(lm1vcov)), grep("splineExp", colnames(lm1vcov)), grep("splineNonExp", colnames(lm1vcov))),
                   c(grep(exposure, colnames(lm1vcov)), grep("splineExp", colnames(lm1vcov)), grep("splineNonExp", colnames(lm1vcov)))]
    v12 <- lm1vcov[grep("splineExp", colnames(lm1vcov)),
                   grep("splineExp", colnames(lm1vcov))]
    v2 <- lm2vcov
    v3 <- (2*sigma(lm2)^4)/lm2$df.residual

    CDEderiv <- c(1, derivPredExp, -derivPredNonExp)
    CDEsigma <- as.matrix(v11)
    CDEv <- c(t(CDEderiv) %*% CDEsigma %*% CDEderiv)
    names(CDEv) <- "CDE (var)"
    NDEderiv <- c(1, derivExpectExp0$dbeta, -derivExpectNonExp0$dbeta,
                  derivExpectExp0$dgamma-derivExpectNonExp0$dgamma,
                  derivExpectExp0$dsigma22-derivExpectNonExp0$dsigma22)
    NDEsigma <- as.matrix(bdiag(v11, v2, v3))
    NDEv <- c(t(NDEderiv) %*% NDEsigma %*% NDEderiv)
    names(NDEv) <- "NDE (var)"
    NIEderiv <- c(derivExpectExp1$dbeta-derivExpectExp0$dbeta,
                  derivExpectExp1$dgamma-derivExpectExp0$dgamma,
                  derivExpectExp1$dsigma22-derivExpectExp0$dsigma22)
    NIEsigma <- as.matrix(bdiag(v12, v2, v3))
    NIEv <- c(t(NIEderiv) %*% NIEsigma %*% NIEderiv)
    names(NIEv) <- "NIE (var)"
  }

  return(list(expMedOutLM = lm1,
              expMedLM = lm2,
              knotSeq = t,
              coefExp = coefExpUpdate,
              coefNonExp = coefNonExpUpdate,
              sigma1 = sigma1,
              residual1 = residual1,
              sigma2 = sigma2,
              residual2 = residual2,
              gamma0 = gamma0,
              gamma1 = gamma1,
              gamma2 = gamma2,
              CDE = CDE,
              NDE = NDE,
              NIE = NIE,
              CDEv = CDEv,
              NDEv = NDEv,
              NIEv = NIEv))
}
