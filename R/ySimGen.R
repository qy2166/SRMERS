#' Generate y using existed X, coefficients of X, variance of random effect and variance of noise
#'
#' @import stats
#' @import dplyr
#' @export ySimGen
#' @param dataset A data frame.
#' @param seed The seed.
#' @param coefVec The coefficient vector of the fixed effect matrix.
#' @param fixMat The fixed effect matrix.
#' @param randomVar The random effect variable in the data set.
#' @param randomSD The SD of random effect.
#' @param noiseSD The SD of noise.
#' @return A data set with generated y.
ySimGen <- function(dataset, seed, coefVec, fixMat, randomVar = NULL, randomSD = NULL, noiseSD){
  if(is.null(randomVar) & is.null(randomSD)){
    # generate noise from a normal distribution
    set.seed(seed)
    dataset$noise <- rnorm(nrow(dataset), 0, noiseSD)

    # generate response vector y
    dataset$ySim <- fixMat %*% coefVec + dataset$noise
    #print("Fixed")
  }
  else if(is.null(randomVar) | is.null(randomSD)){
    stop("Please provide more infomation")
  }
  else{
    # generate random effect from a normal distribution
    randomLen <- dataset %>% group_by_(randomVar) %>% summarise(n = n())
    set.seed(seed)
    randomEff <- rnorm(nrow(randomLen), 0, randomSD)
    randomPrep <- cbind(randomLen, randomEff)
    dataset <- merge(dataset, randomPrep, by = randomVar, all.x = T)

    # generate noise from a normal distribution
    set.seed(seed)
    dataset$noise <- rnorm(nrow(dataset), 0, noiseSD)

    # generate response vector y
    dataset$ySim <- fixMat %*% coefVec + dataset$randomEff + dataset$noise
    #print("Mixed")
  }

  return(dataset)
}
