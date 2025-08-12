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

#' Plot Population
#'
#' Basic plot of the population list
#' @param x Population-list
#' @param type Default "bve" - bv.development, alt: "kinship" - kinship.development(), "pca" - get.pca()
#' @param gen generations to consider
#' @param database groups to consider
#' @param cohorts cohorts to consider
#' @param ... remaining stuff
#' @examples
#' data(ex_pop)
#' plot(ex_pop)
#' @return Summary of the population list including number of individuals, genone length and trait overview
#' @export


plot.population <- function(x, type ="bve",
                               gen = NULL,
                               database = NULL,
                               cohorts = NULL,
                               ...){

  population <- x


  if(length(gen)==0 & length(database)==0 && length(cohorts)==0){
    gen <- 1:length(population$breeding)
  }

  if(type=="bve"){
    bv.development(population, gen = gen, database = database, cohorts = cohorts)
  } else if(type == "kinship"){
    kinship.development(population, gen=gen, database = database, cohorts = cohorts)
  } else if(type == "pca"){
    get.pca(population,  gen = gen, database = database, cohorts = cohorts)
  }

}
