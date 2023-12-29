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

#' Derive class
#'
#' Function to devide the class for each individual
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names (default: FALSE)
#' @param use.all.copy Set to TRUE to extract phenotyping
#' @examples
#' data(ex_pop)
#' get.class(ex_pop, gen=2)
#' @return Class of in gen/database/cohorts selected individuals
#' @export

get.pheno.time <- function(population, database=NULL, gen=NULL, cohorts=NULL, use.id=FALSE, use.all.copy = TRUE){

  database <- get.database(population, gen, database, cohorts)

  n.animals <- sum(database[,4] - database[,3] +1)
  data <- rep(NA, n.animals)
  before <- 0
  names <- numeric(n.animals)

  if(use.all.copy){

    rindex = 1
    for(row in 1:nrow(database)){
      animals <- database[row,]
      if(diff(database[row,3:4])>(-1)){
        for(index in database[row,3]:database[row,4]){
          names[rindex] <- paste(if(animals[2]==1) "M" else "F", index,"_", animals[1],sep="")
          copies <- population$breeding[[animals[1]]][[animals[2]]][[index]][[21]]
          for(rows in 1:nrow(copies)){
            if(sum(population$breeding[[copies[rows,1]]][[copies[rows,2]]][[copies[rows,3]]][[15]]) > 0){
              data[rindex] <- min(data[rindex], population$breeding[[copies[rows,1]]][[42 + copies[rows,2]]][copies[rows,3]], na.rm = TRUE)
            }
          }
          rindex <- rindex + 1
        }
      }
    }


  } else{
    for(row in 1:nrow(database)){
      animals <- database[row,]
      nanimals <- database[row,4] - database[row,3] +1
      if(nanimals>0){
        data[(before+1):(before+nanimals)] <- population$breeding[[animals[1]]][[42+ animals[2]]][animals[3]:animals[4]]
        names[(before+1):(before+nanimals)] <- paste(if(animals[2]==1) "M" else "F", animals[3]:animals[4],"_", animals[1], sep="")
        before <- before + nanimals
      }
    }
  }


  if(use.id){
    names(data) <- get.id(population, database = database)
  } else{
    names(data) <- names
  }

  return(data)
}
