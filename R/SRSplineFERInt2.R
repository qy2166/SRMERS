#' Fit a semi-parametric shape-restricted fixed effect regression spline with factor-by-curve interaction.
#'
#' @import splines2
#' @import stats
#' @import nloptr
#' @export SRSplineFERInt2
#' @param lmFormula A linear model formula.
#' @param dataset A data frame.
#' @param varNameExp The name of the exposure (must be a binary variable).
#' @param varNameMed The name of the mediator (must be a continuous variable).
#' @param knotType The knot type: 1=equal-spaced, 2=quantile, 3=pre-specified.
#' @param preKnot The pre-specified knots.
#' @param nBasis The number of bases.
#' @param Eff0 The spline type for main effect
#' @param Eff1 The spline type for interaction effect
#' @return A list of knots, data set, lm model, summary of lm model, coefficient estimates of bases, covariance of coefficient estimates of bases, projected coefficient estimates, test statistic.
SRSplineFERInt2 <- function(lmFormula, dataset, varNameExp, varNameMed, knotType = 2, preKnot = NULL, nBasis = 5, Eff0 = "I", Eff1 = "I"){
  # step 1: generate knots
  if(knotType == 1){
    knots <- seq(min(dataset[, varNameMed], na.rm = T), max(dataset[, varNameMed], na.rm = T),
                 (max(dataset[, varNameMed], na.rm = T) - min(dataset[, varNameMed], na.rm = T))/(nBasis-1))
    knots <- knots[2:(length(knots)-1)]
  }
  else if(knotType == 2){
    knots <- quantile(dataset[, varNameMed], seq(0, 1, 1/(nBasis-1)), na.rm = T)
    knots <- knots[2:(length(knots)-1)]
  }
  else if(knotType == 3){
    knots <- preKnot
  }
  else{
    stop("please type the correct knot type number")
  }

  # step 2: generate i-spline matrix and c-spline matrix
  iSplineBasisPre <- iSpline(dataset[, varNameMed], knots = knots, degree = 1, intercept = T)
  iSplineBasis0 <- abs(dataset[, varNameExp]-1)*iSplineBasisPre
  colnames(iSplineBasis0) <- paste0("SplineBasis0_", 1:nBasis)
  iSplineBasis1 <- dataset[, varNameExp]*iSplineBasisPre
  colnames(iSplineBasis1) <- paste0("SplineBasis1_", 1:nBasis)

  cSplineBasisPre <- cSpline(dataset[, varNameMed], knots = knots, degree = 1, intercept = T, scale = T)
  cSplineBasisPre <- cbind(dataset[, varNameMed], cSplineBasisPre)
  cSplineBasis0 <- abs(dataset[, varNameExp]-1)*cSplineBasisPre
  colnames(cSplineBasis0) <- paste0("SplineBasis0_", 0:nBasis)
  cSplineBasis1 <- dataset[, varNameExp]*cSplineBasisPre
  colnames(cSplineBasis1) <- paste0("SplineBasis1_", 0:nBasis)

  # step 3: fixed effect regression spline
  if(Eff0 == "I" & Eff1 == "I"){
    datasetNew <- cbind(dataset, iSplineBasis0, iSplineBasis1)
    RSModel <- lm(lmFormula, data = datasetNew)
    RSModelSum <- summary(RSModel)
  }
  else if(Eff0 == "I" & Eff1 == "C"){
    datasetNew <- cbind(dataset, iSplineBasis0, cSplineBasis1)
    RSModel <- lm(lmFormula, data = datasetNew)
    RSModelSum <- summary(RSModel)
  }
  else if(Eff0 == "C" & Eff1 == "I"){
    datasetNew <- cbind(dataset, cSplineBasis0, iSplineBasis1)
    RSModel <- lm(lmFormula, data = datasetNew)
    RSModelSum <- summary(RSModel)
  }
  else if(Eff0 == "C" & Eff1 == "C"){
    datasetNew <- cbind(dataset, cSplineBasis0, cSplineBasis1)
    RSModel <- lm(lmFormula, data = datasetNew)
    RSModelSum <- summary(RSModel)
  }
  else{
    stop("Please enter the correct Eff0 and Eff1")
  }

  # step 4: project paramter onto a cone
  evalFun <- function(theta, thetaHat, sigma){
    objective <- t(theta - thetaHat) %*% solve(sigma) %*% (theta - thetaHat)
    return(objective)
  }
  opts <- list("algorithm" = "NLOPT_LN_BOBYQA",
               "xtol_rel" = 1e-8)
  thetaHat0 <- c(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1") + nBasis - 1)), 1])
  thetaHat1 <- c(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1") + nBasis - 1)), 1])
  sigma0 <- as.matrix(vcov(RSModelSum)[c(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1") + nBasis - 1)),
                                       c(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1") + nBasis - 1))])
  sigma1 <- as.matrix(vcov(RSModelSum)[c(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1") + nBasis - 1)),
                                       c(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1") + nBasis - 1))])

  x0_01 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1") + nBasis - 1)), 1] < 0,
                    RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1") + nBasis - 1)), 1],
                    -0.1)) # negative orthant
  x0_02 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1") + nBasis - 1)), 1] > 0,
                    RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis0_1") + nBasis - 1)), 1],
                    0.1)) # positive orthant

  x0_11 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1") + nBasis - 1)), 1] < 0,
                    RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1") + nBasis - 1)), 1],
                    -0.1)) # negative orthant
  x0_12 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1") + nBasis - 1)), 1] > 0,
                    RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis1_1") + nBasis - 1)), 1],
                    0.1)) # positive orthant

  projResult_0 <- list()
  LRstat_0 <- vector()
  projResult_0[[1]] <- nloptr(x0 = x0_01,
                              eval_f = evalFun,
                              lb = c(rep(-Inf, nBasis)),
                              ub = c(rep(0, nBasis)),
                              opts = opts,
                              thetaHat = thetaHat0,
                              sigma = sigma0)
  LRstat_0[1] <- (t(thetaHat0) %*% solve(sigma0) %*% thetaHat0 - projResult_0[[1]]$objective)/(t(thetaHat0) %*% solve(sigma0) %*% thetaHat0)
  projResult_0[[2]] <- nloptr(x0 = x0_02,
                              eval_f = evalFun,
                              lb = c(rep(0, nBasis)),
                              ub = c(rep(Inf, nBasis)),
                              opts = opts,
                              thetaHat = thetaHat0,
                              sigma = sigma0)
  LRstat_0[2] <- (t(thetaHat0) %*% solve(sigma0) %*% thetaHat0 - projResult_0[[2]]$objective)/(t(thetaHat0) %*% solve(sigma0) %*% thetaHat0)
  names(projResult_0) <- c("neg", "pos")
  names(LRstat_0) <- c("neg", "pos")

  projResult_1 <- list()
  LRstat_1 <- vector()
  projResult_1[[1]] <- nloptr(x0 = x0_11,
                              eval_f = evalFun,
                              lb = c(rep(-Inf, nBasis)),
                              ub = c(rep(0, nBasis)),
                              opts = opts,
                              thetaHat = thetaHat1,
                              sigma = sigma1)
  LRstat_1[1] <- (t(thetaHat1) %*% solve(sigma1) %*% thetaHat1 - projResult_1[[1]]$objective)/(t(thetaHat1) %*% solve(sigma1) %*% thetaHat1)
  projResult_1[[2]] <- nloptr(x0 = x0_12,
                              eval_f = evalFun,
                              lb = c(rep(0, nBasis)),
                              ub = c(rep(Inf, nBasis)),
                              opts = opts,
                              thetaHat = thetaHat1,
                              sigma = sigma1)
  LRstat_1[2] <- (t(thetaHat1) %*% solve(sigma1) %*% thetaHat1 - projResult_1[[2]]$objective)/(t(thetaHat1) %*% solve(sigma1) %*% thetaHat1)
  names(projResult_1) <- c("neg", "pos")
  names(LRstat_1) <- c("neg", "pos")

  return(list(knots = knots,
              dataset = dataset,
              RSModel = RSModel,
              RSModelSum = RSModelSum,
              thetaHat0 = thetaHat0,
              thetaHat1 = thetaHat1,
              sigma0 = sigma0,
              sigma1 = sigma1,
              projResult_0 = projResult_0,
              projResult_1 = projResult_1,
              LRstat_0 = LRstat_0,
              LRstat_1 = LRstat_1))
}
