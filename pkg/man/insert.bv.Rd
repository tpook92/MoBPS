% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insert.bve.R
\name{insert.bv}
\alias{insert.bv}
\title{Manually enter breeding values}
\usage{
insert.bv(
  population,
  bvs,
  na.override = FALSE,
  count = 1,
  count.only.increase = TRUE
)
}
\arguments{
\item{population}{Population list}

\item{bvs}{Matrix of phenotypes to enter (one row per individual with 1 element coding individual name)}

\item{na.override}{Set to TRUE to also enter NA values (Default: FALSE - those entries will be skipped)}

\item{count}{Counting for economic cost calculation (default: 1 - (one observation (for "pheno"), one genotyping (for "bve")))}

\item{count.only.increase}{Set to FALSE to reduce the number of observation for a phenotype to "count" (default: TRUE)}
}
\value{
Population-List with newly entered estimated breeding values
}
\description{
Function to manually enter breeding values
}
\examples{
data(ex_pop)
bv <- get.bv(ex_pop, gen=2)
new.bve <- cbind( colnames(bv), bv[,1]) ## Unrealistic but you do not get better than this!
ex_pop <- insert.bv(ex_pop, bvs=new.bve)
}
