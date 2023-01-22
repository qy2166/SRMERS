#' Calculate the CDE, NDE and NIE (linear models, continuous exposure).
#'
#' @import Matrix
#' @export LRMed2
#' @param data A data frame.
#' @param exposure The name of the exposure (must be a continuous variable).
#' @param mediator The name of the mediator (must be a continuous variable).
#' @param outcome The name of the outcome (must be a continuous variable).
#' @param confounderVec The name vector of the confounders.
#' @param mValue The controlled mediator value for CDE estimation.
#' @param eValueLow The low reference level of exposure.
#' @param eValueHigh The high reference level of exposure.
#' @return A list of exposure-outcome model, exposure-mediator model, CDE, NDE and NIE and their asymptotic variances.
LRMed2 <- function(data, exposure, mediator, outcome, confounderVec, mValue, eValueLow, eValueHigh){
  formula.m <- as.formula(paste0(mediator, "~", exposure, " + ",
                                 paste0(confounderVec, collapse = " + ")))
  formula.y <- as.formula(paste0(outcome, "~", exposure , " + ", mediator, " + ", exposure, ":", mediator, " + ",
                                 paste0(confounderVec, collapse = " + ")))

  model.m <- lm(formula.m, data = data)
  model.y <- lm(formula.y, data = data)

  c = apply(model.m$model[, confounderVec], 2, mean)

  CDE_vander <- (model.y$coefficients[exposure] + model.y$coefficients[paste0(exposure, ":", mediator)]*mValue)*(eValueHigh - eValueLow)
  NDE_vander <- (model.y$coefficients[exposure] +
                   model.y$coefficients[paste0(exposure, ":", mediator)]*(model.m$coefficients %*% c(1, eValueLow, c)))*(eValueHigh - eValueLow)
  NIE_vander <- (model.y$coefficients[mediator]*model.m$coefficients[exposure] +
                   model.y$coefficients[paste0(exposure, ":", mediator)]*model.m$coefficients[exposure]*eValueHigh)*(eValueHigh - eValueLow)

  sigma_vander <- as.matrix(bdiag(vcov(model.m), vcov(model.y)))
  g_CDE_vander <- c(0, 0, rep(0, length(confounderVec)), 0, 1, 0, rep(0, length(confounderVec)), mValue)
  g_NDE_vander <- c(model.y$coefficients[paste0(exposure, ":", mediator)],
                    model.y$coefficients[paste0(exposure, ":", mediator)]*eValueLow,
                    model.y$coefficients[paste0(exposure, ":", mediator)]*c,
                    0, 1, 0, rep(0, length(confounderVec)), model.m$coefficients %*% c(1, eValueLow, c))
  g_NIE_vander <- c(0, model.y$coefficients[mediator] + model.y$coefficients[paste0(exposure, ":", mediator)]*eValueHigh,
                    rep(0, length(confounderVec)), 0, 0, model.m$coefficients[exposure],
                    rep(0, length(confounderVec)), model.m$coefficients[exposure]*eValueHigh)
  CDE_vander_var <- g_CDE_vander %*% sigma_vander %*% g_CDE_vander
  NDE_vander_var <- g_NDE_vander %*% sigma_vander %*% g_NDE_vander
  NIE_vander_var <- g_NIE_vander %*% sigma_vander %*% g_NIE_vander

  CDE <- CDE_vander
  names(CDE) <- "CDE"
  NDE <- NDE_vander[1, 1]
  names(NDE) <- "NDE"
  NIE <- NIE_vander
  names(NIE) <- "NIE"

  CDEv <- CDE_vander_var[1, 1]
  names(CDEv) <- "CDE (var)"
  NDEv <- NDE_vander_var[1, 1]
  names(NDEv) <- "NDE (var)"
  NIEv <- NIE_vander_var[1, 1]
  names(NIEv) <- "NIE (var)"

  return(list(expMedOutLM = model.y,
              expMedLM = model.m,
              CDE = CDE,
              NDE = NDE,
              NIE = NIE,
              CDEv = CDEv,
              NDEv = NDEv,
              NIEv = NIEv))
}
