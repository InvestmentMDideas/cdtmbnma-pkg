\docType{package}
\name{cdtmbnma-package}
\alias{cdtmbnma-package}
\title{Component, Dose, and Time Network Meta-Analysis}
\description{
Bayesian component model-based network meta-analysis for treatment combinations
with explicit component dose-response and dose-dependent interaction surfaces.
The main model supports single-timepoint arm-level networks with continuous or
binary outcomes. A stage-one interface supports a two-component exponential
time-course model for repeated arm means.
}
\details{
The single-timepoint workflow is \code{cdt_data()} followed by \code{cdt_fit()},
or the one-call wrapper \code{cdtmbnma()}. Estimation uses Stan through either
\pkg{cmdstanr} or \pkg{rstan}. Fitted objects have \code{summary()},
\code{coef()}, \code{plot()}, and \code{predict()} methods.
}
\section{Interaction surfaces}{
\describe{
  \item{none}{Additive component effects.}
  \item{bilinear}{One parameter per component pair; the interaction grows with
  the product of normalised component doses.}
  \item{gpdi}{A saturating general pharmacodynamic interaction surface.}
}
}
\author{Tyler Pitre}
\keyword{package}
