% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.parents.R
\name{check.parents}
\alias{check.parents}
\title{Relatedness check between two individuals}
\usage{
check.parents(
  population,
  info.father,
  info.mother,
  max.rel = 2,
  avoid.mating.parent = FALSE,
  still.check = FALSE
)
}
\arguments{
\item{population}{Population list}

\item{info.father}{position of the first parent in the dataset}

\item{info.mother}{position of the second parent in the dataset}

\item{max.rel}{maximal allowed relationship (default: 2, alt: 1 no full-sibs, 0 no half-sibs)}

\item{avoid.mating.parent}{Set to TRUE to avoid matings of an individual to its parents}

\item{still.check}{Internal parameter (avoid.mating.parent check)}
}
\value{
logical with TRUE if relatedness does not excced max.rel / FALSE otherwise.
}
\description{
Internal function to check the relatedness between two individuals
}
\examples{
data(ex_pop)
check.parents(ex_pop, info.father=c(4,1,1,1), info.mother=c(4,2,1,1))
}
