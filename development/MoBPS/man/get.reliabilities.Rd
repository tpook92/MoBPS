% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.reliabilities.R
\name{get.reliabilities}
\alias{get.reliabilities}
\title{Export underlying reliabilities}
\usage{
get.reliabilities(
  population,
  database = NULL,
  gen = NULL,
  cohorts = NULL,
  use.id = FALSE
)
}
\arguments{
\item{population}{Population list}

\item{database}{Groups of individuals to consider for the export}

\item{gen}{Quick-insert for database (vector of all generations to export)}

\item{cohorts}{Quick-insert for database (vector of names of cohorts to export)}

\item{use.id}{Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names (default: FALSE)}
}
\value{
Estimated reliability for BVE for in gen/database/cohorts selected individuals
}
\description{
Function to export underlying reliabilities
}
\examples{
data(ex_pop)
get.reliabilities(ex_pop, gen=2)
}
