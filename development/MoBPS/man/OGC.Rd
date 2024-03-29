% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OGC.R
\name{OGC}
\alias{OGC}
\title{Optimal genetic contribution}
\usage{
OGC(
  A,
  u,
  Q,
  cAc = NA,
  single = TRUE,
  verbose = FALSE,
  max_male = Inf,
  max_female = Inf
)
}
\arguments{
\item{A}{relationship matrix}

\item{u}{breeding values}

\item{Q}{sex indicator}

\item{cAc}{target gain in inbreeding}

\item{single}{If FALSE multiple individuals can be removed at the same type (this is faster but potentially inaccurate!)}

\item{verbose}{Set to FALSE to not display any prints}

\item{max_male}{maximum number of male with positive contributions}

\item{max_female}{maximum number of females with positive contributions}
}
\value{
[[1]] Contributions [[2]] expected inbreeding gain
}
\description{
In this function the OGC selection according to Meuwissen 1997 is performed
}
