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

#' Export parametrization of fixed effects
#'
#' Function to export parametrization of the fixed effects
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names (default: FALSE)
#' @examples
#' data(ex_pop)
#' population <- add.fixed.effects(ex_pop, fixed.effects = cbind(1,5))
#' population <- breeding.diploid(population, heritability = 0.3,
#' fixed.effects.p = rbind(c(1,0), c(0,1)), phenotyping.gen=2)
#' get.fixed.effects.p(population, gen=2)
#' @return Estimated breeding value of in gen/database/cohorts selected individuals
#' @export

get.fixed.effects.p <- function(population, database=NULL, gen=NULL, cohorts=NULL, use.id=TRUE){

  database <- get.database(population, gen, database, cohorts)

  n.animals <- sum(database[,4] - database[,3] +1)
  data <- matrix(0, nrow=n.animals, ncol=ncol(population$info$fixed.effects))
  before <- 1
  names <- numeric(n.animals)
  for(row in 1:nrow(database)){
    animals <- database[row,]
    for(index in database[row,3]:database[row,4]){
      data[before,] <- population$breeding[[animals[1]]][[animals[2]]][[index]][[28]]
      names[before] <- paste(if(animals[2]==1) "M" else "F", index ,"_", animals[1], sep="")
      before <- before +1
    }

  }

  if(use.id){
    rownames(data) <- get.id(population, database = database)
  } else{
    rownames(data) <- names
  }

  return(data)
}
