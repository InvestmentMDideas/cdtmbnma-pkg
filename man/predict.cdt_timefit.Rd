\name{predict.cdt_timefit}
\alias{predict.cdt_timefit}
\title{Predict Long-Term Effects from a Time-Course Fit}
\usage{
\method{predict}{cdt_timefit}(object, newdata, probs = c(0.025, 0.975), ...)
}
\arguments{
\item{object}{A fitted \code{cdt_timefit} object.}
\item{newdata}{Data frame carrying the two component-dose columns used in fitting.}
\item{probs}{Quantiles for credible intervals.}
\item{...}{Unused.}
}
\value{A data frame with doses, posterior mean, posterior standard deviation, and requested quantiles.}
\description{Computes long-term asymptotic relative effects for the stage-one longitudinal model.}
