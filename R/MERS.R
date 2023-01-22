#' P values of for shapes obtained from semi-parametric shape-restricted mixed effects regression splines.
#'
#' @import splines2
#' @import lme4
#' @import lmerTest
#' @import stats
#' @import nloptr
#' @import MASS
#' @export MERS
#' @param y The name of the outcome.
#' @param xMain The name of the main effect.
#' @param xConf The name vector of the confounders.
#' @param xRand The name of the random effect.
#' @param dataset A data frame.
#' @param knotType The knot type: 1=equal-spaced, 2=quantile, 3=pre-specified.
#' @param preKnot The pre-specified knots.
#' @param nBasis The number of bases.
#' @param nIter The number of iterations.
#' @return A list of weights of beta distribution and p-values.
MERS <- function(y, xMain, xConf = NULL, xRand, dataset, knotType = 2, preKnot = NULL, nBasis = 5, nIter){
  lmmformula <- list()
  if(is.null(xConf)){
    lmmformula[[1]] <- as.formula(paste0(y, " ~ ",
                                         paste0("SplineBasis", 1:nBasis, collapse = " + "), " + ",
                                         paste0("(1 | ", xRand, ")")))
    lmmformula[[2]] <- as.formula(paste0(y, " ~ ",
                                         paste0("SplineBasis", 1:nBasis, collapse = " + "), " + ", xMain, " + ",
                                         paste0("(1 | ", xRand, ")")))
  }else{
    lmmformula[[1]] <- as.formula(paste0(y, " ~ ",
                                         paste0("SplineBasis", 1:nBasis, collapse = " + "), " + ",
                                         paste0(xConf, collapse = " + "), " + ",
                                         paste0("(1 | ", xRand, ")")))
    lmmformula[[2]] <- as.formula(paste0(y, " ~ ",
                                         paste0("SplineBasis", 1:nBasis, collapse = "+"), " + ", xMain, " + ",
                                         paste0(xConf, collapse = " + "), " + ",
                                         paste0("(1 | ", xRand, ")")))
  }
  mers <- list()
  mers[[1]] <- iSplineMER(lmmFormula = lmmformula[[1]],
                          dataset = dataset, varName = xMain, knotType = knotType, preKnot = preKnot,
                          nBasis = nBasis, increasing = T)
  mers[[2]] <- iSplineMER(lmmFormula = lmmformula[[1]],
                          dataset = dataset, varName = xMain, knotType = knotType, preKnot = preKnot,
                          nBasis = nBasis, increasing = F)
  mers[[3]] <- cSplineMER(lmmFormula = lmmformula[[2]],
                          dataset = dataset, varName = xMain, knotType = knotType, preKnot = preKnot,
                          nBasis = nBasis, convex = T)
  mers[[4]] <- cSplineMER(lmmFormula = lmmformula[[2]],
                          dataset = dataset, varName = xMain, knotType = knotType, preKnot = preKnot,
                          nBasis = nBasis, convex = F)

  result <- distributionAsymp(nIter = nIter, nBasis = nBasis, sigmai = mers[[1]]$sigma, sigmac = mers[[3]]$sigma,
                              testStatIncr = mers[[1]]$LRstat, testStatDecr = mers[[2]]$LRstat,
                              testStatConv = mers[[3]]$LRstat, testStatConc = mers[[4]]$LRstat)

  return(result)
}
