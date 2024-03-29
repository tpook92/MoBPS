% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.time.point.R
\name{get.time.point}
\alias{get.time.point}
\title{Derive time point}
\usage{
get.time.point(
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
Time point of generation for in gen/database/cohorts selected individuals
}
\description{
Function to devide time point for each individual
}
\examples{
data(ex_pop)
get.time.point(ex_pop, gen=2)
}
