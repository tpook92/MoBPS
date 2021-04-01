'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2017 -- 2020  Torsten Pook

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

#' Estimate effective population size
#'
#' Function to estimate the effective population size
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @examples
#' data(ex_pop)
#' get.effective.size(population=ex_pop, gen=5)
#' @return Estimated effective population size
#' @export


get.effective.size <- function(population, gen=NULL, database=NULL, cohorts=NULL){
  geno <- get.geno(population, gen = gen, database=database, cohorts=cohorts)
  ldinfo <- ld.decay(population, genotype.dataset = geno, type="cm", plot= FALSE)
  effs <- numeric(length(ldinfo[[3]]$x))
  for(index in 1:length(ldinfo[[3]]$x)){
    effs[index] <- effective.size(ldinfo[[3]]$y[index],
                                  dist = ldinfo[[3]]$x[index],
                                  n = ncol(geno))
  }
  return(ceiling(mean(effs[-(1:10)])))
}


#' Estimate effective population size
#'
#' Internal function to estimate the effective population size
#' @param ld ld between markers
#' @param dist distance between markers in Morgan
#' @param n Population size
#'

effective.size <- function(ld, dist, n){

  c <- 0
  for(index in 2*(1:max(5,dist*2))-1){
    c <- c + stats::dpois(index, lambda = dist)
  }

  return( ((1-c)^2 + c^2) / (( ld - 1/n) * (2 * c*(2-c))  ))
}

