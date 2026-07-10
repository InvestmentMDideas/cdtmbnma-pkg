\name{cdt_priors}
\alias{cdt_priors}
\title{Prior Settings for Single-Timepoint Models}
\usage{
cdt_priors(emax_sd = 10, logED50_mean = log(50), logED50_sd = 1,
  int_sd = 5, ref_sd = 10, omega_sd = 2)
}
\arguments{
\item{emax_sd}{Prior standard deviation for component Emax values.}
\item{logED50_mean, logED50_sd}{Prior mean and standard deviation for log half-maximal dose values.}
\item{int_sd}{Prior standard deviation for interaction parameters.}
\item{ref_sd}{Prior standard deviation for per-study reference levels.}
\item{omega_sd}{Half-normal prior scale for the random-effect standard deviation.}
}
\value{A named list of prior hyperparameters.}
\description{Creates prior hyperparameters for \code{cdt_fit()}.}
