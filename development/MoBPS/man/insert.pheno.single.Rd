% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insert.pheno.single.R
\name{insert.pheno.single}
\alias{insert.pheno.single}
\title{Manually enter estimated breeding values}
\usage{
insert.pheno.single(population, phenos, count.only.increase = TRUE)
}
\arguments{
\item{population}{Population list}

\item{phenos}{Matrix of breeding values to enter (one row per individual: 1st column: individual, 2nd column: trait nr, 3rd-nth column records)}

\item{count.only.increase}{Set to FALSE to reduce the number of observation for a phenotype to "count" (default: TRUE)}
}
\value{
Population-List with newly entered estimated breeding values
}
\description{
Function to manually enter estimated breeding values
}
\examples{
data(ex_pop)
bv <- get.bv(ex_pop, gen=2)
new.bve <- cbind( colnames(bv), bv[,1]) ## Unrealistic but you do not get better than this!
ex_pop <- insert.bve(ex_pop, bves=new.bve)
}
