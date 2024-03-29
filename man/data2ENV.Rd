% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data2ENV.R
\name{data2ENV}
\alias{data2ENV}
\title{\emph{data2ENV}: function used to move objects from one environment to another.}
\usage{
data2ENV(data, env)
}
\arguments{
\item{data}{bruSDM data file to be used in the integrated model.}

\item{env}{Environment where the objects should be assigned.}
}
\value{
Assignment of the relevant spatial fields to the specified environment.
}
\description{
Internal function: used to assign objects specified in bruSDM to the dataSDM/blockedCV function environments.
}
