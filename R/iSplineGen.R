#' Generate I-Spline Basis
#'
#' @import splines2
#' @import stats
#' @import graphics
#' @export iSplineGen
#' @param dataset A data frame.
#' @param varName The name of the main effect.
#' @param knotType The knot type: 1=equal-spaced, 2=quantile, 3=pre-specified.
#' @param preKnot The pre-specified knots.
#' @param nBasis The number of bases.
#' @param plot Plot the basis function or not: T=yes, F=no.
#' @return A list of knots, bases, data set.
iSplineGen <- function(dataset, varName, knotType = 1, preKnot = NULL, nBasis = 5, plot = F){
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
  SplineBasisName <- colnames(SplineBasis)
  dataset <- cbind(dataset, SplineBasis)

  if(plot == T){
    SplinePlotData <- data.frame(cbind(dataset[, varName], SplineBasis))
    colnames(SplinePlotData)[1] <- varName
    for(i in 1:ncol(SplineBasis)){
      colnames(SplinePlotData)[i+1] <- paste0("SplineBasis", i)
    }
    SplinePlotData <- SplinePlotData[order(SplinePlotData[, varName]),]
    matplot(SplinePlotData[, varName], SplinePlotData[, 2:(nBasis+1)], type = "l",
            xlab = varName, ylab = "I-spline basis")
    abline(v = knots, lty = 2, col = "gray")
  }
  else if(plot == F){
    print("no plot")
  }
  else{
    stop("plot is a boolean variable")
  }

  return(list(knots = knots,
              SplineBasis = SplineBasis,
              dataset = dataset))
}
