\name{cdt_time_data}
\alias{cdt_time_data}
\title{Assemble a Two-Component Longitudinal Component-Dose Network}
\usage{
cdt_time_data(data, study, arm, time, components, y, se = NULL,
  sd = NULL, n = NULL, ref = NULL, dstar = NULL)
}
\arguments{
\item{data}{A long data frame with one row per study-arm-time observation.}
\item{study}{Name of the study column.}
\item{arm}{Name of the arm column.}
\item{time}{Name of the follow-up time column.}
\item{components}{Character vector naming exactly two component dose columns.}
\item{y, se}{Outcome mean and standard error columns.}
\item{sd, n}{Optional standard deviation and sample-size columns used to derive \code{se}.}
\item{ref}{Optional column flagging the reference arm.}
\item{dstar}{Optional named numeric vector of two dose-normalisation references.}
}
\value{An object of class \code{cdt_time_data}.}
\description{Builds data for the two-component exponential time-course model.}
