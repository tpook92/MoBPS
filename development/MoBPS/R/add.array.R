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
#' @param marker.included Vector with number of SNP entries coding if each marker is on the array (TRUE/FALSE)
#' @param array.name Name of the added array
#' @examples
#' data(ex_pop)
#' population <- add.array(ex_pop, marker.included = c(TRUE, FALSE), array.name="Half-density")
#' @return Population list
#' @export

add.array <- function(population, marker.included = TRUE,
                      array.name = NULL){

  if(length(array.name)==0){
    array.name = paste0("Array_", length(population$info$array.name)+1)
  }
  if(length(marker.included)<sum(population$info$snp)){
    marker.included <- rep(marker.included, length.out = sum(population$info$snp))
  }

  population$info$array.name = c(population$info$array.name, array.name)
  population$info$array.markers[[length(population$info$array.markers)+1]] = marker.included
  population$info$array.is_subset =  c(population$info$array.is_subset , prod(marker.included)!=1)

  return(population)
}




