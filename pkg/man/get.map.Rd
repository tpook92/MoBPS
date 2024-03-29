% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.map.R
\name{get.map}
\alias{get.map}
\title{Map generation}
\usage{
get.map(population, use.snp.nr = FALSE, morgan.position.per.chromosome = TRUE)
}
\arguments{
\item{population}{Population list}

\item{use.snp.nr}{Set to TRUE to display SNP number and not SNP name}

\item{morgan.position.per.chromosome}{Set to FALSE to Morgan position continuously over the genome}
}
\value{
Genomic map of the population list
}
\description{
Function to derive the genomic map for a given population list
}
\examples{
data(ex_pop)
map <- get.map(ex_pop)
}
