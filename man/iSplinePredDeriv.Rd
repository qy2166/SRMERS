% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iSplinePredDeriv.R
\name{iSplinePredDeriv}
\alias{iSplinePredDeriv}
\title{Calculate the derivatives of coefficients of a fixed effect regression spline (monotonicity).}
\usage{
iSplinePredDeriv(t, theta, m)
}
\arguments{
\item{t}{The knot sequence vector in exposure-outcome model.}

\item{theta}{The coefficient vector of I-spline bases in exposure-outcome model.}

\item{m}{The mediator value within min(t) and max(t).}
}
\value{
The derivatives of beta.
}
\description{
Calculate the derivatives of coefficients of a fixed effect regression spline (monotonicity).
}
