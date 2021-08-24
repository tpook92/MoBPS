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

#' Set defaults
#'
#' Set default parameter values in breeding.diploid
#' @param population Population list
#' @param parameter.name Number of traits (If more than traits via real.bv.X use traits with no directly underlying QTL)
#' @param parameter.value Genetic variance of traits with no underlying QTL
#' @param parameter.remove Remove a specific previously generated parameter default
#' @param reset.all Set to TRUE to remove all prior parameter values
#' @return Population-list with one or more additional new traits
#' @examples
#' data(ex_pop)
#' population <- set.default(ex_pop, parameter.name="heritability", parameter.value=0.3)
#' @export


set.default <- function(population,
                        parameter.name=NULL,
                        parameter.value=NULL,
                        parameter.remove=NULL,
                        reset.all=FALSE){


  if(reset.all){
    population$info$default.parameter.name = NULL
    population$info$default.parameter.value = list()
  }

  if(length(parameter.name)>0){

    replace <- which(population$info$default.parameter.name == parameter.name)

    if(length(replace)==1){

      population$info$default.parameter.value[[replace]] <- parameter.value
    } else{

      population$info$default.parameter.name <- c(population$info$default.parameter.name, parameter.name)

      population$info$default.parameter.value[[ length(population$info$default.parameter.value)+1]] <- parameter.value
    }


  }


  if(length(parameter.remove)==1){

    replace <- which(population$info$default.parameter == parameter.name)

    if(length(replace)==1){

      population$info$default.parameter.name <- population$info$default.parameter.name[-replace]
      population$info$default.parameter.value[[replace]] <- NULL
    } else{

      warning("nothing to remove here!")

    }
  } else if(length(parameter.remove)>1){
    warning("you cant remove multiple parameter simulatanously with parameter.remove")
  }




  return(population)
}
