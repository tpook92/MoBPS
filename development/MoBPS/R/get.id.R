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

#' Derive ID on an individual
#'
#' Function to derive the internal ID given to each individual
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names
#' @param keep.order To not change order of individuals when ids are provided (default: FALSE)
#' @examples
#' data(ex_pop)
#' get.id(ex_pop, gen=2)
#' @return Individual ID for in gen/database/cohorts selected individuals
#' @export

get.id <- function(population, database=NULL, gen=NULL, cohorts=NULL, use.id=FALSE, keep.order = FALSE){

  database <- get.database(population, gen, database, cohorts, keep.order = keep.order)

  n.animals <- sum(database[,4] - database[,3] +1)
  data <- rep(0, n.animals)
  before <- 0
  names <- numeric(n.animals)
  for(row in 1:nrow(database)){
    animals <- database[row,]
    nanimals <- database[row,4] - database[row,3] +1
    if(nanimals>0){
      data[(before+1):(before+nanimals)] <- population$breeding[[animals[1]]][[14+ animals[2]]][animals[3]:animals[4]]
      names[(before+1):(before+nanimals)] <- paste(if(animals[2]==1) "M" else "F", animals[3]:animals[4],"_", animals[1], sep="")
      before <- before + nanimals
    }
  }
  if(use.id){
    names(data) <- get.id(population, database = database)
  } else{
    names(data) <- names
  }

  return(data)
}
