\name{predict.cdtmbnma}
\alias{predict.cdtmbnma}
\title{Predict Relative Effects at Component-Dose Combinations}
\usage{
\method{predict}{cdtmbnma}(object, newdata, probs = c(0.025, 0.975), ...)
}
\arguments{
\item{object}{A fitted \code{cdtmbnma} object.}
\item{newdata}{Data frame carrying the component dose columns used in fitting.}
\item{probs}{Quantiles for credible intervals.}
\item{...}{Unused.}
}
\value{A data frame with doses, posterior mean, posterior standard deviation, and requested quantiles.}
\description{
Computes posterior relative effects against the all-zero reference at arbitrary
component-dose combinations.
}
