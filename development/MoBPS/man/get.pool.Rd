% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.pool.R
\name{get.pool}
\alias{get.pool}
\title{Export underlying true breeding values}
\usage{
get.pool(
  population,
  gen = NULL,
  database = NULL,
  cohorts = NULL,
  ids = NULL,
  plot = FALSE,
  import.position.calculation = NULL,
  decodeOriginsU = decodeOriginsR
)
}
\arguments{
\item{population}{Population list}

\item{gen}{Quick-insert for database (vector of all generations to export)}

\item{database}{Groups of individuals to consider for the export}

\item{cohorts}{Quick-insert for database (vector of names of cohorts to export)}

\item{ids}{AA}

\item{plot}{Set TRUE to generate a visualization of genetic origins}

\item{import.position.calculation}{Function to calculate recombination point into adjacent/following SNP}

\item{decodeOriginsU}{Used function for the decoding of genetic origins [[5]]/[[6]]}
}
\value{
Genomic value of in gen/database/cohorts selected individuals
}
\description{
Function to export underlying true breeding values
}
\examples{
data(ex_pop)
get.bv(ex_pop, gen=2)
}
