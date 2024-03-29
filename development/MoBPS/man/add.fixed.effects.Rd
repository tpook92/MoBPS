% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.fixed.effects.R
\name{add.fixed.effects}
\alias{add.fixed.effects}
\title{Add a trait as a linear combination of other traits}
\usage{
add.fixed.effects(population, fixed.effects, replace = FALSE)
}
\arguments{
\item{population}{population list}

\item{fixed.effects}{Matrix containing fixed effects (p x k -matrix with p being the number of traits and k being number of fixed effects; default: not fixed effects (NULL))}

\item{replace}{Set to TRUE to delete previously added fixed effects}
}
\value{
Population list

Population list
}
\description{
Function to create an additional trait that is the results of a linear combination of the other traits
}
\examples{
data(ex_pop)
population <- add.fixed.effects(ex_pop, fixed.effects = matrix(c(3,5), nrow=1))
}
