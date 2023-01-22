#' P values of for shapes obtained from semi-parametric shape-restricted fixed effect regression splines with factor-by-curve interaction.
#'
#' @import splines2
#' @import stats
#' @import nloptr
#' @import MASS
#' @export FERSInt
#' @param y The name of the outcome.
#' @param xExp The name of the exposure (must be a binary variable).
#' @param xMed The name of the mediator (must be a continuous variable).
#' @param xConf The name vector of the confounders.
#' @param dataset A data frame.
#' @param knotType The knot type: 1=equal-spaced, 2=quantile, 3=pre-specified.
#' @param preKnot The pre-specified knots.
#' @param nBasis The number of bases.
#' @param nIter The number of iterations.
#' @return A list of weights of beta distribution and p-values for both exposure groups.
FERSInt <- function(y, xExp, xMed, xConf = NULL, dataset, knotType = 2, preKnot = NULL, nBasis = 5, nIter){
  lmformula <- list()
  if(is.null(xConf)){
    lmformula[[1]] <- as.formula(paste0(y, " ~ ",
                                        paste0("SplineBasis_", 1:nBasis, collapse = " + "), " + ",
                                        paste0("SplineBasisInt_", 1:nBasis, collapse = " + ")))
    lmformula[[2]] <- as.formula(paste0(y, " ~ ",
                                        paste0("SplineBasis_", 1:nBasis, collapse = " + "), " + ",
                                        paste0("SplineBasisInt_", 1:nBasis, collapse = " + "), " + ",
                                        "SplineBasis_0 + SplineBasisInt_0"))
  }else{
    lmformula[[1]] <- as.formula(paste0(y, " ~ ",
                                        paste0("SplineBasis_", 1:nBasis, collapse = " + "), " + ",
                                        paste0("SplineBasisInt_", 1:nBasis, collapse = " + "), " + ",
                                        paste0(xConf, collapse = " + ")))
    lmformula[[2]] <- as.formula(paste0(y, " ~ ",
                                        paste0("SplineBasis_", 1:nBasis, collapse = "+"), " + ",
                                        paste0("SplineBasisInt_", 1:nBasis, collapse = " + "), " + ",
                                        "SplineBasis_0 + SplineBasisInt_0", " + ",
                                        paste0(xConf, collapse = " + ")))
  }
  fers <- list()
  fers[[1]] <- SRSplineFERInt(lmformula[[1]],
                              dataset = dataset,
                              varNameExp = xExp, varNameMed = xMed, knotType = 2, preKnot = preKnot, nBasis = nBasis,
                              EffMain = "I", EffInt = "I")
  fers[[2]] <- SRSplineFERInt(lmformula[[2]],
                              dataset = dataset,
                              varNameExp = xExp, varNameMed = xMed, knotType = 2, preKnot = preKnot, nBasis = nBasis,
                              EffMain = "C", EffInt = "C")

  resultMain <- distributionAsymp(nIter = nIter, nBasis = nBasis, sigmai = fers[[1]]$sigmaMain, sigmac = fers[[2]]$sigmaMain,
                                  testStatIncr = fers[[1]]$LRstat_Main[2], testStatDecr = fers[[1]]$LRstat_Main[1],
                                  testStatConv = fers[[2]]$LRstat_Main[2], testStatConc = fers[[2]]$LRstat_Main[1])

  resultAdd <- distributionAsymp(nIter = nIter, nBasis = nBasis, sigmai = fers[[1]]$sigmaAdd, sigmac = fers[[2]]$sigmaAdd,
                                 testStatIncr = fers[[1]]$LRstat_Add[2], testStatDecr = fers[[1]]$LRstat_Add[1],
                                 testStatConv = fers[[2]]$LRstat_Add[2], testStatConc = fers[[2]]$LRstat_Add[1])

  return(list(resultMain, resultAdd))
}
