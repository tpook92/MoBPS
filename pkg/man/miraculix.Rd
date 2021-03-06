% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miraculix.R
\name{miraculix}
\alias{miraculix}
\title{Add miraculix-coding for genotypes}
\usage{
miraculix(population)
}
\arguments{
\item{population}{Population list}
}
\value{
Population list
}
\description{
Internal function to store genotypes bit-wise
}
\examples{
# This is only relevant with the package miraculix is installed and used
population <- creating.diploid(nsnp=100, nindi=50, miraculix=FALSE)
\donttest{population <- miraculix(population)}
}
