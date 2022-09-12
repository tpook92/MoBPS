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
#' @param fixed.effects Matrix containing fixed effects (p x k -matrix with p being the number of traits and k being number of fixed effects; default: not fixed effects (NULL))
#' @param replace Set to TRUE to delete previously added fixed effects
#' @return Population list
#' @examples
#' data(ex_pop)
#' population <- add.fixed.effects(ex_pop, fixed.effects = matrix(c(3,5), nrow=1))
#' @return Population list
#' @export
#'
add.fixed.effects <- function(population, fixed.effects, replace = FALSE){

  if(is.matrix(fixed.effects)){
    if(length(fixed.effects)< population$info$bv.nr){
      fixed.effects <- rep(fixed.effects, length.out = population$info$bv.nr)
    }
    fixed.effects <- matrix(fixed.effects, nrow= population$info$bv.nr, ncol = length(fixed.effects) / population$info$bv.nr)
  }

  if(replace){
    population$info$fixed.effects <- fixed.effects
  } else{
    population$info$fixed.effects <- cbind(population$info$fixed.effects, fixed.effects)
  }



  temp1 <- c(rep(0, ncol(population$info$fixed.effects)))
  for(gen in 1:length(population$breeding)){
    for(sex in 1:2){
      if(length(population$breeding[[gen]][[sex]])>0){
        for(index in 1:length(population$breeding[[gen]][[sex]])){
          if(replace){
            population$breeding[[gen]][[sex]][[index]][[28]] <- temp1
          } else{
            population$breeding[[gen]][[sex]][[index]][[28]] <- c(population$breeding[[gen]][[sex]][[index]][[28]], rep(0, ncol(population$info$fixed.effects) - length(population$breeding[[gen]][[sex]][[index]][[28]])))
          }
        }
      }
    }
  }


  return(population)
}




