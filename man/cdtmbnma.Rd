\name{cdtmbnma}
\alias{cdtmbnma}
\title{One-Call Component Dose-Response Network Meta-Analysis}
\usage{
cdtmbnma(data, study, components, outcome = c("continuous", "binary"),
  y = NULL, se = NULL, sd = NULL, n = NULL,
  events = NULL, n_binary = NULL, ref = NULL,
  dstar = NULL, interactions = NULL,
  interaction = c("bilinear", "none", "gpdi"),
  priors = cdt_priors(), newdata = NULL, backend = "auto", chains = 4,
  iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
  seed = 1, refresh = 0, ...)
}
\arguments{
\item{data, study, components, outcome, y, se, sd, n, events, n_binary, ref, dstar, interactions}{Passed to \code{cdt_data()}.}
\item{interaction, priors, newdata, backend, chains, iter_warmup, iter_sampling, adapt_delta, seed, refresh, ...}{Passed to \code{cdt_fit()}.}
}
\value{An object of class \code{cdtmbnma}.}
\description{
Convenience wrapper that builds the design with \code{cdt_data()} and fits it
with \code{cdt_fit()}.
}
\examples{
\donttest{
# Fitting needs a Stan backend, and the Stan model is compiled on first use.
if (requireNamespace("rstan", quietly = TRUE)) {
  fit <- cdtmbnma(sacval_example(), study = "study",
                  components = c("d_sac", "d_val"), outcome = "continuous",
                  y = "y", se = "se", backend = "rstan",
                  chains = 2, iter_warmup = 500, iter_sampling = 500)
  summary(fit)
}
}
}
