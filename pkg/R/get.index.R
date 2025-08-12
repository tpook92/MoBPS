'#
  Authors
Torsten Pook, torsten.pook@wur.nl

Copyright (C) 2017 -- 2025  Torsten Pook

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

#' Putting together indices for GxE / multi trait
#'
#' Function to put together indices for GxE / multi trait
#' @param population Population list
#' @param traits Which traits to include in the index (all weight with factor 1)
#' @param locations Which locations to include in the index (all weight weight factor 1)
#' @param trait.weights Vector with a weight per trait
#' @param location.weights Vector weight a weight per location
#' @return Index
#' @examples
#' population = creating.diploid(nsnp =1000, nindi = 100)
#' population = creating.trait(population, n.additive = c(10,10), n.location=3, replace.traits = TRUE)
#' get.index(population, trait.weights = c(1,2), location.weights = c(1,2,3))
#' @export

get.index <- function(population, traits = NULL, locations = NULL, trait.weights = NULL, location.weights = NULL){

  if(length(trait.weights)==0 & length(traits)==0){
    trait.weights = rep(1, max(population$info$trait.nr))
  }

  if(length(traits)>0){
    trait.weights = rep(0, max(population$info$trait.nr))
    trait.weights[traits] = 1
  }

  if(length(location.weights)==0 & length(locations)==0){
    location.weights = rep(1, max(population$info$trait.location))
  }

  if(length(locations)>0){
    location.weights = rep(0, max(population$info$trait.location))
    location.weights[locations] = 1
  }

  index = location.weights[population$info$trait.location] * trait.weights[population$info$trait.nr]

  return(index)
}
