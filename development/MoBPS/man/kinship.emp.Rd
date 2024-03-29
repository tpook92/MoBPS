% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinship.emp.R
\name{kinship.emp}
\alias{kinship.emp}
\title{Empirical kinship}
\usage{
kinship.emp(
  population = NULL,
  animals = NULL,
  gen = NULL,
  database = NULL,
  cohorts = NULL,
  sym = FALSE,
  use.id = FALSE
)
}
\arguments{
\item{population}{Population list}

\item{animals}{List of animals to compute kinship for}

\item{gen}{Quick-insert for database (vector of all generations to export)}

\item{database}{Groups of individuals to consider for the export}

\item{cohorts}{Quick-insert for database (vector of names of cohorts to export)}

\item{sym}{If True derive matrix entries below principle-diagonal}

\item{use.id}{Set TRUE to use animal IDs for column/row-names in the output matrix}
}
\value{
Empirical kinship matrix (IBD-based since Founders)
}
\description{
Function to compute empirical kinship for a set of individuals)
}
\examples{
data(ex_pop)
kinship <- kinship.emp(population=ex_pop, database=cbind(2,1,1,25))
}
