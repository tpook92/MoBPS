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

#' Exclude individuals from a database
#'
#' Function to exclude individuals from a database
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param remove.database Groups of individuals to remove from the database (same IDs!)
#' @param remove.gen Generations of individuals to remove from the database (same IDs!)
#' @param remove.cohorts Cohorts of individuals to remove from the database (same IDs!)
#' @return Database excluding removals
#' @examples
#' data(ex_pop)
#' database <- group.diff(ex_pop, gen=1, remove.database=cbind(1,1))
#' @export

group.diff <- function(population, database=NULL, gen=NULL, cohorts=NULL, remove.gen = NULL, remove.database = NULL, remove.cohorts= NULL){


  all_id <- get.id(population, gen=gen, database=database, cohorts=cohorts)

  old_database <- get.database(population, gen, database, cohorts)

  if(length(remove.gen)==0 && length(remove.database)==0 && length(remove.cohorts)==0){
    return(old_database)
  }

  remove_id <- unique(sort(get.id(population, gen=remove.gen, database=remove.database, cohorts=remove.cohorts)))


  new_database <- matrix(0, ncol=4, nrow=length(all_id))

  nr1 <- 1
  nr2 <- 1
  for(index in 1:nrow(old_database)){
    activ <- old_database[index,]
    for(index2 in activ[3]:activ[4]){
      if(sum(remove_id==all_id[nr1])==0){
        new_database[nr2,] <- c(activ[1], activ[2], index2,index2)
        nr2 <- nr2 +1
      }
      nr1 <- nr1 +1
    }
  }

  if(nr2==1){
    new_database = matrix(nrow = 0, ncol = 4)
  } else{
    new_database <- new_database[1:(nr2-1),]
  }


  new_database <- get.database(population, database = new_database)

  return(new_database)

}
