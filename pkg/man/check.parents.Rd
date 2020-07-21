% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.parents.R
\name{check.parents}
\alias{check.parents}
\title{Derive pedigree}
\usage{
check.parents(population, info.father, info.mother, max.rel = 2)
}
\arguments{
\item{population}{Population list}

\item{info.father}{position of the first parent in the dataset}

\item{info.mother}{position of the second parent in the dataset}

\item{max.rel}{maximal allowed relationship (default: 2, alt: 1 no full-sibs, 0 no half-sibs)}
}
\value{
Pedigree-file for in gen/database/cohorts selected individuals
}
\description{
Derive pedigree for selected individuals
}
\examples{
data(ex_pop)
get.pedigree(ex_pop, gen=2)
}