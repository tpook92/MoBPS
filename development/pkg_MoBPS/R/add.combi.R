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

#' Add a trait as a linear combination of other traits
#'
#' Function to create an additional trait that is the results of a linear combination of the other traits
#' @param population population list
#' @param trait trait nr. for which to implement a combination of other traits
#' @param combi.weights Weights (only linear combinations of other traits are allowed!)
#' @param trait.name  Name of the trait generated
#' @return Population list
#' @examples
#' data(ex_pop)
#' population <- creating.trait(ex_pop, n.additive = 100)
#' population <- add.combi(population, trait = 3, combi.weights = c(1,5))
#' @return Population list
#' @export
#'
add.combi <- function(population, trait, combi.weights, trait.name = NULL){

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
  if(length(trait.name)==1){
    population$info$trait.name[trait] <- trait.name
  }

  combi.weights <- c(combi.weights, rep(0, population$info$bv.nr - length(combi.weights)))

  if(combi.weights[trait]!=0){
    stop("No weighting on trait itself allowed!")
  }

  population$info$combi.weights[[trait]] <- combi.weights

  population$info$bv.calculated <- FALSE
  population <- breeding.diploid(population)

  return(population)
}




