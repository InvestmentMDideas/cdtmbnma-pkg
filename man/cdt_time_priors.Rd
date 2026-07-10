\name{cdt_time_priors}
\alias{cdt_time_priors}
\title{Prior Settings for the Two-Component Time-Course Model}
\usage{
cdt_time_priors(ref_E_sd = 5, ref_lograte_mean = -3,
  ref_lograte_sd = 1, emax_sd = 10, logED50_mean = 0,
  logED50_sd = 1, eta_sd = 5, rate_sd = 1,
  omega_E_sd = 2, omega_k_sd = 0.5)
}
\arguments{
\item{ref_E_sd}{Prior standard deviation for study-specific reference asymptotes.}
\item{ref_lograte_mean, ref_lograte_sd}{Prior mean and standard deviation for study-specific reference log rates.}
\item{emax_sd}{Prior standard deviation for component Emax values.}
\item{logED50_mean, logED50_sd}{Prior mean and standard deviation for log half-maximal dose values.}
\item{eta_sd}{Prior standard deviation for the bilinear interaction.}
\item{rate_sd}{Prior standard deviation for component log-rate dose effects.}
\item{omega_E_sd, omega_k_sd}{Half-normal prior scales for random-effect standard deviations.}
}
\value{A named list of prior hyperparameters.}
\description{Creates prior hyperparameters for \code{cdt_time_fit()}.}
