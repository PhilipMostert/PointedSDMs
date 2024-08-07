% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataSet.R
\name{dataSet}
\alias{dataSet}
\title{Internal function used to standardize datasets, as well as assign metadata.}
\usage{
dataSet(
  datapoints,
  datanames,
  coords = c("CoordLoc1", "CoordLoc2"),
  proj,
  pointcovnames,
  paresp,
  countsresp,
  trialname,
  speciesname,
  marks,
  marktrialname,
  markfamily,
  temporalvar,
  offsetname
)
}
\arguments{
\item{datapoints}{A list of datasets as sf objects}

\item{datanames}{A vector of the names of the datasets.}

\item{coords}{Names of the coordinates used in the model.}

\item{proj}{The projection reference system used in the model.}

\item{pointcovnames}{Name of the point covariates used in the model.}

\item{paresp}{Name of the response variable used by the presence absence datasets.}

\item{countsresp}{Name of the response variable used by the counts data.}

\item{trialname}{Name of the trial variable used by the presence absence datasets.}

\item{speciesname}{Name of the species name variable.}

\item{marks}{Name of the marks considered in the model.}

\item{marktrialname}{Name of the trial variable used by the binomial marks.}

\item{markfamily}{A vector describing the distribution of the marks.}

\item{temporalvar}{Name of the temporal variable.}

\item{offsetname}{Name of the offset column in the datasets.}
}
\value{
A list of relevant metadata
}
\description{
Internal function used to assist in structuring the data.
}
