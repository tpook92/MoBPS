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

#' Check for copies in genotype calculation
#'
#' Check for copies in genotype calculation
#' @param population population-list
#' @param to.compute database to check for copied animals
#' @return database with first copies
#' @export

check.copy = function(population, to.compute){

  if(length(to.compute)==3){

    activ = to.compute
    first_copy = population$breeding[[activ[1]]][[activ[2]]][[activ[3]]][[21]][1:3]

    if(prod(activ==first_copy)!=1){
      to.compute = matrix(first_copy, ncol = 3)
    }

  } else{
    for(index in 1:nrow(to.compute)){

      activ = to.compute[index,]
      first_copy = population$breeding[[activ[1]]][[activ[2]]][[activ[3]]][[21]][1:3]

      if(prod(activ==first_copy)!=1){
        to.compute[index,] = first_copy
      }
    }
  }


  return(to.compute)

}
