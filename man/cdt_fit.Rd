\name{cdt_fit}
\alias{cdt_fit}
\title{Fit a Component Dose-Response Network Meta-Analysis}
\usage{
cdt_fit(data, interaction = c("bilinear", "none", "gpdi"),
  priors = cdt_priors(), newdata = NULL, backend = "auto", chains = 4,
  iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
  seed = 1, refresh = 0, ...)
}
\arguments{
\item{data}{A \code{cdt_data} object.}
\item{interaction}{Interaction surface: \code{"bilinear"}, \code{"none"}, or \code{"gpdi"}.}
\item{priors}{A named list from \code{cdt_priors()}.}
\item{newdata}{Optional data frame of component-dose combinations for Stan generated-quantities predictions.}
\item{backend}{One of \code{"auto"}, \code{"cmdstanr"}, or \code{"rstan"}.}
\item{chains, iter_warmup, iter_sampling}{Sampler settings.}
\item{adapt_delta}{Target acceptance probability passed to the Stan backend.}
\item{seed}{Random seed.}
\item{refresh}{Console refresh interval.}
\item{...}{Passed to the backend sampler.}
}
\value{An object of class \code{cdtmbnma}.}
\description{
Compiles the Stan model on first use and samples the posterior.
}
\details{
The additive model uses \code{interaction = "none"}. The bilinear model uses one
interaction parameter per co-occurring component pair. The \code{"gpdi"} model
uses a saturating surface adapted from the general pharmacodynamic interaction
model.
}
\examples{
sv <- sacval_example()
d <- cdt_data(sv, "study", c("d_sac", "d_val"), outcome = "continuous",
              y = "y", se = "se", dstar = c(d_sac = 200, d_val = 320))

\donttest{
# Fitting needs a Stan backend, and the Stan model is compiled on first use.
if (requireNamespace("rstan", quietly = TRUE)) {
  fit <- cdt_fit(d, interaction = "bilinear", backend = "rstan",
                 chains = 2, iter_warmup = 500, iter_sampling = 500)
  summary(fit)
}
}
}
