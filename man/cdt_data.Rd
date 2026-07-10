\name{cdt_data}
\alias{cdt_data}
\title{Assemble a Component Dose-Response Network}
\usage{
cdt_data(data, study, components, outcome = c("continuous", "binary"),
  y = NULL, se = NULL, sd = NULL, n = NULL,
  events = NULL, n_binary = NULL, ref = NULL,
  dstar = NULL, interactions = NULL)
}
\arguments{
\item{data}{A data frame with one row per study arm.}
\item{study}{Name of the column identifying the study.}
\item{components}{Character vector of component dose column names. A value of zero marks the component absent.}
\item{outcome}{Either \code{"continuous"} or \code{"binary"}.}
\item{y}{For a continuous outcome, the arm mean column.}
\item{se}{For a continuous outcome, the arm standard error column. Supply this or both \code{sd} and \code{n}.}
\item{sd, n}{For a continuous outcome, columns of arm standard deviation and sample size.}
\item{events, n_binary}{For a binary outcome, event-count and sample-size columns.}
\item{ref}{Optional column flagging exactly one within-study reference arm per study.}
\item{dstar}{Optional named numeric vector of dose-normalisation references for the bilinear surface.}
\item{interactions}{\code{NULL}, or a list of length-2 character vectors naming component pairs to model.}
}
\value{An object of class \code{cdt_data}.}
\description{
Turns an arm-level data frame into the design consumed by \code{cdt_fit()}.
Every component named in \code{components} must have a dose column, and zero
means the component is absent from that arm.
}
\examples{
df <- data.frame(study = c("s1", "s1", "s1"), A = c(0, 1, 1),
                 B = c(0, 0, 1), y = c(0, -1, -2), se = c(1, 1, 1))
d <- cdt_data(df, study = "study", components = c("A", "B"),
              outcome = "continuous", y = "y", se = "se")
d
}
