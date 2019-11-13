'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2017 -- 2018  Torsten Pook

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
'#

#' BV standardization
#'
#' Function to get mean and genetic variance of a trait to a fixed value
#' @param population Population list
#' @param mean.target Target mean
#' @param var.target Target variance
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @export


bv.standardization <- function(population, mean.target=100, var.target=10, gen=NULL, database=NULL, cohorts=NULL){

  n_traits <- population$info$bv.nr
  if(length(mean.target)<n_traits) mean.target <- rep(mean.target, length.out = n_traits)
  if(length(var.target)<n_traits) var.target <- rep(var.target, length.out = n_traits)
  if(length(gen)==0 && length(database)==0 && length(cohorts)==0){
    gen <- nrow(population$info$size)
  }
  database <- get.database(population, gen, database, cohorts)
  ## Variance Standardization
  for(index in 1:n_traits){
    new_var <- var.target[index]

    if(population$info$bv.calculated==FALSE){
      population <- breeding.diploid(population, verbose=FALSE)
    }

    var_test <- stats::var(get.bv(population, database= database)[index,])
    test1 <- TRUE
    if(length(population$info$real.bv.add[[index]])>0){
      population$info$real.bv.add[[index]][,3:5] <- population$info$real.bv.add[[index]][,3:5] * sqrt(  new_var / var_test)
      test1 <- FALSE
    }
    if(length(population$info$real.bv.mult[[index]])>0){
      population$info$real.bv.mult[[index]][,5:13] <- population$info$real.bv.mult[[index]][,5:13] * sqrt(  new_var / var_test)
      test1 <- FALSE
    }
    if(test1 && verbose) cat("You entered a trait without quantitative loci. Is this intentional?\n")

  }
  population$info$bv.calculated <- FALSE

  ## Mean Standardization
  for(index in 1:n_traits){

    if(population$info$bv.calculated==FALSE){
      population <- breeding.diploid(population, verbose=FALSE)
    }

    mean_test <- mean(get.bv(population, database = database)[index,])


    population$info$base.bv[index] <- mean.target[index] + population$info$base.bv[index] - mean_test

  }
  population$info$bv.calculated <- FALSE

  return(population)
}
