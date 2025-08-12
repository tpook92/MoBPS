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

#' Number of generations
#'
#' Function to calculate the number of generations in the population list
#' @param population Population list
#' @return Numeric value
#' @examples
#' data(ex_pop)
#' get.ngen(ex_pop)
#' @export


get.ngen <- function(population){
  return(length(population$breeding))
}

#' Number of individuals in each generation
#'
#' Function to calculate the number of individuals per generation
#' @param population Population list
#' @return matrix with numeric values
#' @examples
#' data(ex_pop)
#' get.size(ex_pop)
#' @export


get.size <- function(population){
  return(population$info$size)
}

#' Number of traits
#'
#' Function to calculate the number of traits in the population list
#' @param population Population list
#' @return Numeric value
#' @examples
#' data(ex_pop)
#' get.ngen(ex_pop)
#' @export


get.ntrait <- function(population){
  return(length(population$info$trait.name))
}

#' Name of traits
#'
#' Function to export trait names in the population list
#' @param population Population list
#' @return Numeric value
#' @examples
#' data(ex_pop)
#' get.ngen(ex_pop)
#' @export


get.trait.name <- function(population){
  return((population$info$trait.name))
}
