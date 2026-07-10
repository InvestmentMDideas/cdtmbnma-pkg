\name{cdt_time_fit}
\alias{cdt_time_fit}
\title{Fit the Two-Component Component-Dose-Time Model}
\usage{
cdt_time_fit(data, priors = cdt_time_priors(), newdata = NULL,
  backend = "auto", chains = 4, iter_warmup = 1000,
  iter_sampling = 1000, adapt_delta = 0.95, seed = 1,
  refresh = 0, ...)
}
\arguments{
\item{data}{A \code{cdt_time_data} object.}
\item{priors}{A named list from \code{cdt_time_priors()}.}
\item{newdata}{Optional data frame of two-component dose combinations.}
\item{backend}{One of \code{"auto"}, \code{"cmdstanr"}, or \code{"rstan"}.}
\item{chains, iter_warmup, iter_sampling}{Sampler settings.}
\item{adapt_delta}{Target acceptance probability.}
\item{seed}{Random seed.}
\item{refresh}{Console refresh interval.}
\item{...}{Passed to the backend sampler.}
}
\value{An object of class \code{cdt_timefit}.}
\description{Fits the stage-one longitudinal Stan model with an exponential time-course.}
