#' Fit a semi-parametric shape-restricted fixed effect regression spline (monotonicity).
#'
#' @import splines2
#' @import stats
#' @import nloptr
#' @export iSplineFER
#' @param lmFormula A linear model formula.
#' @param dataset A data frame.
#' @param varName The name of the main effect.
#' @param knotType The knot type: 1=equal-spaced, 2=quantile, 3=pre-specified.
#' @param preKnot The pre-specified knots.
#' @param nBasis The number of bases.
#' @param increasing Increasing shape or Decreasing shape: T=Increasing, F=Decreasing
#' @return A list of knots, data set, lm model, summary of lm model, coefficient estimates of bases, covariance of coefficient estimates of bases, combined coefficient estimates, test statistic.
iSplineFER <- function(lmFormula, dataset, varName, knotType = 2, preKnot = NULL, nBasis = 5, increasing = T){
  # step 1: generate knots
  if(knotType == 1){
    knots <- seq(min(dataset[, varName], na.rm = T), max(dataset[, varName], na.rm = T),
                 (max(dataset[, varName], na.rm = T) - min(dataset[, varName], na.rm = T))/(nBasis-1))
    knots <- knots[2:(length(knots)-1)]
  }
  else if(knotType == 2){
    knots <- quantile(dataset[, varName], seq(0, 1, 1/(nBasis-1)), na.rm = T)
    knots <- knots[2:(length(knots)-1)]
  }
  else if(knotType == 3){
    knots <- preKnot
  }
  else{
    stop("please type the correct knot type number")
  }

  # step 2: generate i-spline matrix
  SplineBasis <- iSpline(dataset[, varName], knots = knots, degree = 1, intercept = T)
  for(i in 1:ncol(SplineBasis)){
    colnames(SplineBasis)[i] <- paste0("SplineBasis", i)
  }
  datasetNew <- cbind(dataset, SplineBasis)

  # step 3: fixed effect regression spline
  FERSModel <- lm(formula = lmFormula, data = datasetNew)
  FERSModelSum <- summary(FERSModel)

  # step 4: project paramter onto a cone
  if(increasing == T){
    evalFun <- function(theta, thetaHat, sigma){
      objective <- t(theta - thetaHat) %*% solve(sigma) %*% (theta - thetaHat)
      return(objective)
    }
    x0 <- ifelse(FERSModelSum$coefficients[c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1)), 1] > 0,
                 FERSModelSum$coefficients[c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1)), 1],
                 0.1)
    opts <- list("algorithm" = "NLOPT_LN_BOBYQA",
                 "xtol_rel" = 1e-8)
    thetaHat <- FERSModelSum$coefficients[c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1)), 1]
    sigma <- as.matrix(vcov(FERSModelSum)[c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1)),
                                          c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1))])
    projResult <- nloptr(x0 = x0,
                         eval_f = evalFun,
                         lb = c(rep(0, nBasis)),
                         ub = c(rep(Inf, nBasis)),
                         opts = opts,
                         thetaHat = thetaHat,
                         sigma = sigma)
    # LRstat <- t(thetaHat) %*% solve(sigma) %*% thetaHat - projResult$objective
    LRstat <- (t(thetaHat) %*% solve(sigma) %*% thetaHat - projResult$objective)/(t(thetaHat) %*% solve(sigma) %*% thetaHat)
    # LRstat <- (t(thetaHat) %*% solve(sigma) %*% thetaHat - projResult$objective)/(t(residuals(FERSModel)) %*% residuals(FERSModel) + t(thetaHat) %*% solve(sigma) %*% thetaHat)
  }
  else if(increasing == F){
    evalFun <- function(theta, thetaHat, sigma){
      objective <- t(theta - thetaHat) %*% solve(sigma) %*% (theta - thetaHat)
      return(objective)
    }
    x0 <- ifelse(FERSModelSum$coefficients[c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1)), 1] < 0,
                 FERSModelSum$coefficients[c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1)), 1],
                 -0.1)
    opts <- list("algorithm" = "NLOPT_LN_BOBYQA",
                 "xtol_rel" = 1e-8)
    thetaHat <- FERSModelSum$coefficients[c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1)), 1]
    sigma <- as.matrix(vcov(FERSModelSum)[c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1)),
                                          c(which(rownames(FERSModelSum$coefficients) == "SplineBasis1"):(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis - 1))])
    projResult <- nloptr(x0 = x0,
                         eval_f = evalFun,
                         lb = c(rep(-Inf, nBasis)),
                         ub = c(rep(0, nBasis)),
                         opts = opts,
                         thetaHat = thetaHat,
                         sigma = sigma)
    # LRstat <- t(thetaHat) %*% solve(sigma) %*% thetaHat - projResult$objective
    LRstat <- (t(thetaHat) %*% solve(sigma) %*% thetaHat - projResult$objective)/(t(thetaHat) %*% solve(sigma) %*% thetaHat)
    # LRstat <- (t(thetaHat) %*% solve(sigma) %*% thetaHat - projResult$objective)/(t(residuals(FERSModel)) %*% residuals(FERSModel) + t(thetaHat) %*% solve(sigma) %*% thetaHat)
  }
  else{
    stop("increasing is a boolean variable")
  }

  # step 5: output the parameter
  if((which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis) <= nrow(FERSModelSum$coefficients)){
    iSplineFERPar <- c(FERSModelSum$coefficients[c(1:(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") - 1)), 1],
                       projResult$solution,
                       FERSModelSum$coefficients[c((which(rownames(FERSModelSum$coefficients) == "SplineBasis1") + nBasis):nrow(FERSModelSum$coefficients)), 1])
  }
  else{
    iSplineFERPar <- c(FERSModelSum$coefficients[c(1:(which(rownames(FERSModelSum$coefficients) == "SplineBasis1") - 1)), 1],
                       projResult$solution)
  }
  names(iSplineFERPar)[2:(2+nBasis-1)] <- paste0("SplineBasis", 1:nBasis)
  names(iSplineFERPar)[1] <- "(Intercept)"

  return(list(knots = knots,
              dataset = dataset,
              FERSModel = FERSModel,
              FERSModelSum = FERSModelSum,
              thetaHat = thetaHat,
              sigma = sigma,
              projResult = projResult,
              SplineFERPar = iSplineFERPar,
              LRstat = LRstat))
}
