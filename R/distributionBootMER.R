#' Bootstrap a null distribution of the test statistic for mixed effect model
#'
#' @import splines2
#' @import lme4
#' @import lmerTest
#' @import stats
#' @import nloptr
#' @export distributionBootMER
#' @param nBoot The number of bootstrap iterations.
#' @param seed The seed.
#' @param parametric The method of blup and noise bootstrap: T=resample from estimated blup and noise, F=resample from normal
#' @param MERSObject1 The lmm model under null distribution.
#' @param MERSObject2 The lmm model fitted using iSplineMER or cSplineMER
#' @param randomFormula The formula used to generate random effect matrix (in the form of ~random+0)
#' @param splineRule iSplineMER or cSplineMER
#' @param lmmFormula A linear mixed effect model formula (same as lmmFormula in iSplineMER or cSplineFMER).
#' @param varName The name of the main effect (same as varName in iSplineMER or cSplineMER).
#' @param knotType The knot type: 1=equal-spaced and 2=quantile (same as knotType in iSplineMER or cSplineMER).
#' @param nBasis The number of bases (same as nBasis in iSplineMER or cSplineMER).
#' @param increasing Increasing shape or Decreasing shape: T=Increasing, F=Decreasing
#' @param convex Convex shape or Concave shape: T=Convex, F=Concave
#' @param testStat The test statistic generated from iSplineMER or cSplineMER
#' @return A list of null distribution of test statistic and p value.
distributionBootMER <- function(nBoot, seed, parametric = T, MERSObject1, MERSObject2, randomFormula, splineRule, lmmFormula, varName, knotType = 1, nBasis = 5, increasing = T, convex = T, testStat = NULL){
  # fixed effect generation
  fixCoef <- summary(MERSObject1)$coefficients[, 1]
  if(ncol(MERSObject1@frame) > 1){
    fixMat <- cbind(intercept = 1, MERSObject1@frame[, 2:nrow(summary(MERSObject1)$coefficients)])
  }
  else{
    fixMat <- rep(1, nrow(MERSObject1@frame))
  }
  fixEff <- as.matrix(fixMat) %*% fixCoef

  randomMat <- model.matrix(randomFormula, data = MERSObject1@frame)

  LRstat <- vector()

  set.seed(seed)
  for(i in 1:nBoot){
    if(parametric == T){
      blup <- ranef(MERSObject2$MERSModel)
      randomCoef <- as.matrix(sample(blup[[1]][, 1], replace = T))
      randomEff <- randomMat %*% randomCoef

      # noise
      residual <- residuals(MERSObject2$MERSModel)
      noise <- as.matrix(sample(residual, replace = T))
    }
    else{
      # random effect generation
      randomCoef <- rnorm(ncol(randomMat), 0, sqrt(MERSObject2$MERSModelSum$varcor$cluster[1, 1]))
      randomEff <- randomMat %*% randomCoef

      # noise
      noise <- rnorm(nrow(MERSObject2$MERSModel@frame), 0, MERSObject2$MERSModelSum$sigma)
    }

    dataNew <- MERSObject2$dataset
    dataNew$yBoot <- fixEff + randomEff + noise

    LRstat[i] <- splineRule(lmmFormula = lmmFormula,
                            dataset = dataNew, varName = varName, knotType = knotType, nBasis = nBasis, increasing = increasing, convex = convex)$LRstat
    print(paste0("running the ", i, "th iteration"))
  }

  if(!is.null(testStat)){
    p_value <- sum(testStat < LRstat)/nBoot
  }else{
    p_value <- NULL
  }

  return(list(nullDist = LRstat,
              pValue = p_value))
}
