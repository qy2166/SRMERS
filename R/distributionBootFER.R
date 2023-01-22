#' Bootstrap a null distribution of the test statistic for fixed effect model
#'
#' @import splines2
#' @import stats
#' @import nloptr
#' @export distributionBootFER
#' @param nBoot The number of bootstrap iterations.
#' @param seed The seed.
#' @param parametric The method of noise bootstrap: T=resample from estimated noise, F=resample from normal
#' @param FERSObject1 The lm model under null distribution.
#' @param FERSObject2 The lm model fitted using iSplineFER or cSplineFER
#' @param splineRule iSplineFER or cSplineFER
#' @param lmFormula A linear model formula (same as lmFormula in iSplineFER or cSplineFER).
#' @param varName The name of the main effect (same as varName in iSplineFER or cSplineFER).
#' @param knotType The knot type: 1=equal-spaced and 2=quantile (same as knotType in iSplineFER or cSplineFER).
#' @param nBasis The number of bases (same as nBasis in iSplineFER or cSplineFER).
#' @param increasing Increasing shape or Decreasing shape: T=Increasing, F=Decreasing
#' @param convex Convex shape or Concave shape: T=Convex, F=Concave
#' @param testStat The test statistic generated from iSplineFER or cSplineFER
#' @return A list of null distribution of test statistic and p value.
distributionBootFER <- function(nBoot, seed, parametric = T, FERSObject1, FERSObject2, splineRule, lmFormula, varName, knotType = 1, nBasis = 5, increasing = T, convex = T, testStat = NULL){
  # fixed effect generation
  fixCoef <- summary(FERSObject1)$coefficients[, 1]
  if(ncol(FERSObject1$model) > 1){
    fixMat <- cbind(intercept = 1, FERSObject1$model[, 2:nrow(summary(FERSObject1)$coefficients)])
  }
  else{
    fixMat <- rep(1, nrow(FERSObject1$model))
  }
  fixEff <- as.matrix(fixMat) %*% fixCoef

  LRstat <- vector()

  set.seed(seed)
  for(i in 1:nBoot){
    if(parametric == T){
      # noise
      residual <- residuals(FERSObject2$FERSModel)
      noise <- as.matrix(sample(residual, replace = T))
    }
    else{
      # noise
      noise <- rnorm(nrow(FERSObject2$FERSModel$model), 0, FERSObject2$FERSModelSum$sigma)
    }

    dataNew <- FERSObject2$dataset
    dataNew$yBoot <- fixEff + noise

    LRstat[i] <- splineRule(lmFormula = lmFormula,
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
