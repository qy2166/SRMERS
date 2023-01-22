#' Fit a semi-parametric shape-restricted fixed effect regression spline with factor-by-curve interaction.
#'
#' @import splines2
#' @import stats
#' @import nloptr
#' @export SRSplineFERInt
#' @param lmFormula A linear model formula.
#' @param dataset A data frame.
#' @param varNameExp The name of the exposure (must be a binary variable).
#' @param varNameMed The name of the mediator (must be a continuous variable).
#' @param knotType The knot type: 1=equal-spaced, 2=quantile, 3=pre-specified.
#' @param preKnot The pre-specified knots.
#' @param nBasis The number of bases.
#' @param EffMain The spline type for main effect
#' @param EffInt The spline type for interaction effect
#' @return A list of knots, data set, lm model, summary of lm model, coefficient estimates of bases, covariance of coefficient estimates of bases, projected coefficient estimates, test statistic.
SRSplineFERInt <- function(lmFormula, dataset, varNameExp, varNameMed, knotType = 2, preKnot = NULL, nBasis = 5, EffMain = "I", EffInt = "I"){
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
  iSplineBasis <- iSplineBasisPre
  colnames(iSplineBasis) <- paste0("SplineBasis_", 1:nBasis)
  iSplineBasisInt <- dataset[, varNameExp]*iSplineBasisPre
  colnames(iSplineBasisInt) <- paste0("SplineBasisInt_", 1:nBasis)

  cSplineBasisPre <- cSpline(dataset[, varNameMed], knots = knots, degree = 1, intercept = T, scale = T)
  cSplineBasisPre <- cbind(dataset[, varNameMed], cSplineBasisPre)
  cSplineBasis <- cSplineBasisPre
  colnames(cSplineBasis) <- paste0("SplineBasis_", 0:nBasis)
  cSplineBasisInt <- dataset[, varNameExp]*cSplineBasisPre
  colnames(cSplineBasisInt) <- paste0("SplineBasisInt_", 0:nBasis)

  # step 3: fixed effect regression spline
  if(EffMain == "I" & EffInt == "I"){
    datasetNew <- cbind(dataset, iSplineBasis, iSplineBasisInt)
    RSModel <- lm(lmFormula, data = datasetNew)
    RSModelSum <- summary(RSModel)
  }
  else if(EffMain == "I" & EffInt == "C"){
    datasetNew <- cbind(dataset, iSplineBasis, cSplineBasisInt)
    RSModel <- lm(lmFormula, data = datasetNew)
    RSModelSum <- summary(RSModel)
  }
  else if(EffMain == "C" & EffInt == "I"){
    datasetNew <- cbind(dataset, cSplineBasis, iSplineBasisInt)
    RSModel <- lm(lmFormula, data = datasetNew)
    RSModelSum <- summary(RSModel)
  }
  else if(EffMain == "C" & EffInt == "C"){
    datasetNew <- cbind(dataset, cSplineBasis, cSplineBasisInt)
    RSModel <- lm(lmFormula, data = datasetNew)
    RSModelSum <- summary(RSModel)
  }
  else{
    stop("Please enter the correct EffMain and EffInt")
  }

  # step 4: project paramter onto a cone
  evalFun <- function(theta, thetaHat, sigma){
    objective <- t(theta - thetaHat) %*% solve(sigma) %*% (theta - thetaHat)
    return(objective)
  }
  opts <- list("algorithm" = "NLOPT_LN_BOBYQA",
               "xtol_rel" = 1e-8)
  thetaHatMain <- c(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)), 1])
  thetaHatInt <- c(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)), 1])
  sigmaMain <- as.matrix(vcov(RSModelSum)[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)),
                                          c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1))])
  sigmaInt <- as.matrix(vcov(RSModelSum)[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)),
                                         c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1))])

  x0_Main1 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)), 1] < 0,
                       RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)), 1],
                       -0.1)) # negative orthant
  x0_Main2 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)), 1] > 0,
                       RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)), 1],
                       0.1)) # positive orthant

  x0_Int1 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)), 1] < 0,
                      RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)), 1],
                      -0.1)) # negative orthant
  x0_Int2 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)), 1] > 0,
                      RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)), 1],
                      0.1)) # positive orthant

  projResult_Main <- list()
  LRstat_Main <- vector()
  projResult_Main[[1]] <- nloptr(x0 = x0_Main1,
                                 eval_f = evalFun,
                                 lb = c(rep(-Inf, nBasis)),
                                 ub = c(rep(0, nBasis)),
                                 opts = opts,
                                 thetaHat = thetaHatMain,
                                 sigma = sigmaMain)
  LRstat_Main[1] <- (t(thetaHatMain) %*% solve(sigmaMain) %*% thetaHatMain - projResult_Main[[1]]$objective)/(t(thetaHatMain) %*% solve(sigmaMain) %*% thetaHatMain)
  projResult_Main[[2]] <- nloptr(x0 = x0_Main2,
                                 eval_f = evalFun,
                                 lb = c(rep(0, nBasis)),
                                 ub = c(rep(Inf, nBasis)),
                                 opts = opts,
                                 thetaHat = thetaHatMain,
                                 sigma = sigmaMain)
  LRstat_Main[2] <- (t(thetaHatMain) %*% solve(sigmaMain) %*% thetaHatMain - projResult_Main[[2]]$objective)/(t(thetaHatMain) %*% solve(sigmaMain) %*% thetaHatMain)
  names(projResult_Main) <- c("neg", "pos")
  names(LRstat_Main) <- c("neg", "pos")

  projResult_Int <- list()
  LRstat_Int <- vector()
  projResult_Int[[1]] <- nloptr(x0 = x0_Int1,
                                eval_f = evalFun,
                                lb = c(rep(-Inf, nBasis)),
                                ub = c(rep(0, nBasis)),
                                opts = opts,
                                thetaHat = thetaHatInt,
                                sigma = sigmaInt)
  LRstat_Int[1] <- (t(thetaHatInt) %*% solve(sigmaInt) %*% thetaHatInt - projResult_Int[[1]]$objective)/(t(thetaHatInt) %*% solve(sigmaInt) %*% thetaHatInt)
  projResult_Int[[2]] <- nloptr(x0 = x0_Int2,
                                eval_f = evalFun,
                                lb = c(rep(0, nBasis)),
                                ub = c(rep(Inf, nBasis)),
                                opts = opts,
                                thetaHat = thetaHatInt,
                                sigma = sigmaInt)
  LRstat_Int[2] <- (t(thetaHatInt) %*% solve(sigmaInt) %*% thetaHatInt - projResult_Int[[2]]$objective)/(t(thetaHatInt) %*% solve(sigmaInt) %*% thetaHatInt)
  names(projResult_Int) <- c("neg", "pos")
  names(LRstat_Int) <- c("neg", "pos")

  projResult_Add <- list()
  LRstat_Add <- vector()
  thetaHatAdd <- NA
  sigmaAdd <- NA
  if((EffMain == "I" & EffInt == "I")|(EffMain == "C" & EffInt == "C")){
    thetaHatAdd <- thetaHatMain + thetaHatInt
    sigmaSpline <- as.matrix(vcov(RSModelSum)[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1),
                                                which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)),
                                              c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1),
                                                which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1))])
    sigmaMult <- cbind(diag(1, nrow = nBasis, ncol = nBasis), diag(1, nrow = nBasis, ncol = nBasis))
    sigmaAdd <- sigmaMult %*% sigmaSpline %*% t(sigmaMult)

    x0_Add1 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)), 1] +
                          RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)), 1] < 0,
                        RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)), 1] +
                          RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)), 1],
                        -0.1)) # negative orthant
    x0_Add2 <- c(ifelse(RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)), 1] +
                          RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)), 1] > 0,
                        RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasis_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasis_1") + nBasis - 1)), 1] +
                          RSModelSum$coefficients[c(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1"):(which(rownames(RSModelSum$coefficients) == "SplineBasisInt_1") + nBasis - 1)), 1],
                        0.1)) # positive orthant

    projResult_Add[[1]] <- nloptr(x0 = x0_Add1,
                                  eval_f = evalFun,
                                  lb = c(rep(-Inf, nBasis)),
                                  ub = c(rep(0, nBasis)),
                                  opts = opts,
                                  thetaHat = thetaHatAdd,
                                  sigma = sigmaAdd)
    LRstat_Add[1] <- (t(thetaHatAdd) %*% solve(sigmaAdd) %*% thetaHatAdd - projResult_Add[[1]]$objective)/(t(thetaHatAdd) %*% solve(sigmaAdd) %*% thetaHatAdd)
    projResult_Add[[2]] <- nloptr(x0 = x0_Add2,
                                  eval_f = evalFun,
                                  lb = c(rep(0, nBasis)),
                                  ub = c(rep(Inf, nBasis)),
                                  opts = opts,
                                  thetaHat = thetaHatAdd,
                                  sigma = sigmaAdd)
    LRstat_Add[2] <- (t(thetaHatAdd) %*% solve(sigmaAdd) %*% thetaHatAdd - projResult_Add[[2]]$objective)/(t(thetaHatAdd) %*% solve(sigmaAdd) %*% thetaHatAdd)
    names(projResult_Add) <- c("neg", "pos")
    names(LRstat_Add) <- c("neg", "pos")
  }

  return(list(knots = knots,
              dataset = dataset,
              RSModel = RSModel,
              RSModelSum = RSModelSum,
              thetaHatMain = thetaHatMain,
              thetaHatInt = thetaHatInt,
              thetaHatAdd = thetaHatAdd,
              sigmaMain = sigmaMain,
              sigmaInt = sigmaInt,
              sigmaAdd = sigmaAdd,
              projResult_Main = projResult_Main,
              projResult_Int = projResult_Int,
              projResult_Add = projResult_Add,
              LRstat_Main = LRstat_Main,
              LRstat_Int = LRstat_Int,
              LRstat_Add = LRstat_Add))
}
