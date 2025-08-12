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

#' Relatedness check between two individuals
#'
#' Internal function to check the relatedness between two individuals
#' @param population Population list
#' @param info.father position of the first parent in the dataset
#' @param info.mother position of the second parent in the dataset
#' @param max.rel maximal allowed relationship (default: 2, alt: 1 no full-sibs, 0 no half-sibs)
#' @param avoid.mating.parent Set to TRUE to avoid matings of an individual to its parents
#' @param still.check Internal parameter (avoid.mating.parent check)
#' @examples
#' data(ex_pop)
#' check.parents(ex_pop, info.father=c(4,1,1,1), info.mother=c(4,2,1,1))
#' @return logical with TRUE if relatedness does not excced max.rel / FALSE otherwise.
#' @export


check.parents <- function(population, info.father, info.mother, max.rel=2,
                          avoid.mating.parent = FALSE, still.check = FALSE){

  if(max.rel==2 && avoid.mating.parent){

    if(still.check){
      p1 <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[7]]
      p2 <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[8]]
      p3 <- population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[7]]
      p4 <- population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[8]]

      check <- prod(population$breeding[[p1[1]]][[p1[2]]][[p1[3]]][[21]][1,] == population$breeding[[p3[1]]][[p3[2]]][[p3[3]]][[21]][1,]) +
        prod(population$breeding[[p2[1]]][[p2[2]]][[p2[3]]][[21]][1,] == population$breeding[[p4[1]]][[p4[2]]][[p4[3]]][[21]][1,])

      return(c(TRUE, check))
    } else{
      return(c(TRUE, 0))
    }

  }  else{
    p1 <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[7]]
    p2 <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[8]]
    p3 <- population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[7]]
    p4 <- population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[8]]

    check <- prod(population$breeding[[p1[1]]][[p1[2]]][[p1[3]]][[21]][1,] == population$breeding[[p3[1]]][[p3[2]]][[p3[3]]][[21]][1,]) +
      prod(population$breeding[[p2[1]]][[p2[2]]][[p2[3]]][[21]][1,] == population$breeding[[p4[1]]][[p4[2]]][[p4[3]]][[21]][1,])

    check1 = check2 = 0
    if(avoid.mating.parent){

      check1 <- prod(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[21]][1,] == population$breeding[[p3[1]]][[p3[2]]][[p3[3]]][[21]][1,]) +
        prod(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[21]][1,] == population$breeding[[p4[1]]][[p4[2]]][[p4[3]]][[21]][1,])

      check2 <- prod(population$breeding[[p1[1]]][[p1[2]]][[p1[3]]][[21]][1,] == population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[21]][1,]) +
        prod(population$breeding[[p2[1]]][[p2[2]]][[p2[3]]][[21]][1,] == population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[21]][1,])

    }

    valid = (check<=max.rel) && check1==0 && check2 == 0

    return(c(valid, check))
  }

}
