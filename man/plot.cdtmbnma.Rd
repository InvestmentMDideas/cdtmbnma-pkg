\name{plot.cdtmbnma}
\alias{plot.cdtmbnma}
\title{Plot Component Dose-Response Curves}
\usage{
\method{plot}{cdtmbnma}(x, ngrid = 60, probs = c(0.025, 0.975), ...)
}
\arguments{
\item{x}{A fitted \code{cdtmbnma} object.}
\item{ngrid}{Number of dose points per curve.}
\item{probs}{Two quantiles for the credible band.}
\item{...}{Passed to \code{plot()}.}
}
\value{Invisibly, a list of plotted grids.}
\description{Draws marginal component dose-response curves with posterior credible bands.}
