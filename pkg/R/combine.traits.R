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

#' Combine traits
#'
#' Function to combine traits in the BVE
#' @param population Population list
#' @param combine.traits Vector containing the traits (numbers) to combine into a joined trait
#' @param combine.name Name of the combined trait
#' @param remove.combine Remove a selected previously generated combined trait
#' @param remove.all Set TRUE to remove all previously generated combined traits
#' @return Population-list
#' @examples
#' population <- creating.diploid(nsnp=100, nindi=100, n.additive = c(50,50))
#' population <- combine.traits(population, combine.traits=1:2)
#' population <- breeding.diploid(population, bve=TRUE, phenotyping.gen=1, heritability=0.3)
#' @export



combine.traits <- function(population, combine.traits=NULL, combine.name=NULL, remove.combine=NULL, remove.all=FALSE){

  if(remove.all){
    population$info$trait.combine.name = NULL
    population$info$trait.combine.included = list()
  }

  if(length(combine.traits)>0){

    population$info$trait.combine.included[[length(population$info$trait.combine.included)+1]] <- combine.traits

    if(length(combine.name)!=1){
      combine.name <- paste0("Combined Trait ", length(population$info$trait.combine.included))
    }

    population$info$trait.combine.name <- c(population$info$trait.combine.name, combine.name)

  }

  in_combine <- numeric(population$info$bv.nr)

  for(index in 1:length(population$info$trait.combine.included)){
    in_combine[population$info$trait.combine.included[[index]]] <-     in_combine[population$info$trait.combine.included[[index]]] + 1
  }

  if(sum(in_combine>1)>0){
    stop("A Trait can only be part of one combined trait!")
  }

  return(population)

}
