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

#' Add a genotyping array
#'
#' Function to add a genotyping array for the population
#' @param population population list
#' @param trait trait nr. for which to implement a combination of other traits
#' @param combi.weights Weights (only linear combinations of other traits are allowed!)
#' @examples
#' Hier sollte ihre Werbung stehen!!!
#' @return Population list
#' @export

add.combi <- function(population, trait, combi.weights){

  if(trait > population$info$bv.nr){
    if(trait == (population$info$bv.nr+1)){
      population <- creating.trait(population, bv.total=trait, verbose=FALSE)
    } else{
      stop("Illegal trait to enter a trait combination")
    }
  }
  population$info$is.combi[trait] <- TRUE
  population$info$is.maternal[trait] <- FALSE
  population$info$is.paternal[trait] <- FALSE

  combi.weights <- c(combi.weights, rep(0, population$info$bv.nr - length(combi.weights)))

  if(combi.weights[trait]!=0){
    stop("No weighting on trait itself allowed!")
  }

  population$info$combi.weights[[trait]] <- combi.weights


  return(population)
}




