#' P values of for shapes obtained from semi-parametric shape-restricted fixed effect regression splines.
#'
#' @import splines2
#' @import stats
#' @import nloptr
#' @import MASS
#' @export FERS
#' @param y The name of the outcome.
#' @param xMain The name of the main effect.
#' @param xConf The name vector of the confounders.
#' @param dataset A data frame.
#' @param knotType The knot type: 1=equal-spaced, 2=quantile, 3=pre-specified.
#' @param preKnot The pre-specified knots.
#' @param nBasis The number of bases.
#' @param nIter The number of iterations.
#' @return A list of weights of beta distribution and p-values.
FERS <- function(y, xMain, xConf = NULL, dataset, knotType = 2, preKnot = NULL, nBasis = 5, nIter){
  lmformula <- list()
  if(is.null(xConf)){
    lmformula[[1]] <- as.formula(paste0(y, " ~ ", 
                                        paste0("SplineBasis", 1:nBasis, collapse = " + ")))
    lmformula[[2]] <- as.formula(paste0(y, " ~ ", 
                                        paste0("SplineBasis", 1:nBasis, collapse = " + "), " + ", xMain))
  }else{
    lmformula[[1]] <- as.formula(paste0(y, " ~ ", 
                                        paste0("SplineBasis", 1:nBasis, collapse = " + "), " + ", 
                                        paste0(xConf, collapse = " + ")))
    lmformula[[2]] <- as.formula(paste0(y, " ~ ", 
                                        paste0("SplineBasis", 1:nBasis, collapse = "+"), " + ", xMain, " + ", 
                                        paste0(xConf, collapse = " + ")))
  }
  fers <- list()
  fers[[1]] <- iSplineFER(lmFormula = lmformula[[1]],
                          dataset = dataset, varName = xMain, knotType = knotType, preKnot = preKnot,
                          nBasis = nBasis, increasing = T)
  fers[[2]] <- iSplineFER(lmFormula = lmformula[[1]],
                          dataset = dataset, varName = xMain, knotType = knotType, preKnot = preKnot,
                          nBasis = nBasis, increasing = F)
  fers[[3]] <- cSplineFER(lmFormula = lmformula[[2]],
                          dataset = dataset, varName = xMain, knotType = knotType, preKnot = preKnot,
                          nBasis = nBasis, convex = T)
  fers[[4]] <- cSplineFER(lmFormula = lmformula[[2]],
                          dataset = dataset, varName = xMain, knotType = knotType, preKnot = preKnot,
                          nBasis = nBasis, convex = F)
  
  result <- distributionAsymp(nIter = nIter, nBasis = nBasis, sigmai = fers[[1]]$sigma, sigmac = fers[[3]]$sigma,
                              testStatIncr = fers[[1]]$LRstat, testStatDecr = fers[[2]]$LRstat,
                              testStatConv = fers[[3]]$LRstat, testStatConc = fers[[4]]$LRstat)
  
  return(result)
}