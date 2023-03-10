% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iSplineMER.R
\name{iSplineMER}
\alias{iSplineMER}
\title{Fit a semi-parametric shape-restricted mixed effects regression spline (monotonicity).}
\usage{
iSplineMER(
  lmmFormula,
  dataset,
  varName,
  knotType = 2,
  preKnot = NULL,
  nBasis = 5,
  increasing = T
)
}
\arguments{
\item{lmmFormula}{A linear mixed effect model formula.}

\item{dataset}{A data frame.}

\item{varName}{The name of the main effect.}

\item{knotType}{The knot type: 1=equal-spaced, 2=quantile, 3=pre-specified.}

\item{preKnot}{The pre-specified knots.}

\item{nBasis}{The number of bases.}

\item{increasing}{Increasing shape or Decreasing shape: T=Increasing, F=Decreasing}
}
\value{
A list of knots, data set, lmm model, summary of lmm model, coefficient estimates of bases, covariance of coefficient estimates of bases, combined coefficient estimates, test statistic.
}
\description{
Fit a semi-parametric shape-restricted mixed effects regression spline (monotonicity).
}
