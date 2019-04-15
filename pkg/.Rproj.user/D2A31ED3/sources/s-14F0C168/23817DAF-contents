'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2017 -- 2018  Torsten Pook

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

#' Derive genetic origins
#'
#' Function to derive genetic origin
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @export

get.recombi <- function(population, database=NULL, gen=NULL, cohorts=NULL){

  database <- get.database(population, gen, database, cohorts)

  n.animals <- sum(database[,4] - database[,3] +1)

  data <- list()
  rindex <- 1

  for(row in 1:nrow(database)){
    animals <- database[row,]
    for(index in database[row,3]:database[row,4]){
      data[[rindex]] <-  list()

      colnamed <- paste(if(animals[2]==1) "M" else "W", index,"_", animals[1],sep="")
      data[[rindex]][[1]] <- population$breeding[[animals[1]]][[animals[2]]][[index]][[1]]
      data[[rindex]][[2]] <- population$breeding[[animals[1]]][[animals[2]]][[index]][[2]]
      if(population$info$miraculix){
        if (requireNamespace("miraculix", quietly = TRUE)) {
          decode <- miraculix::decodeOrigins
        } else{
          decode <- decodeOriginsR
        }
      } else{
        decode <- decodeOriginsR
      }
      if(population$info$miraculix){
        male_ori <- matrix(0, nrow=length(population$breeding[[animals[1]]][[animals[2]]][[index]][[5]]), ncol=4)
        female_ori <- matrix(0, nrow=length(population$breeding[[animals[1]]][[animals[2]]][[index]][[6]]), ncol=4)
        for(index2 in 1:nrow(male_ori)){
          male_ori[index2,] <- decode(population$breeding[[animals[1]]][[animals[2]]][[index]][[5]],index2)
        }
        for(index2 in 1:nrow(female_ori)){
          female_ori[index2,] <- decode(population$breeding[[animals[1]]][[animals[2]]][[index]][[6]],index2)
        }
        data[[rindex]][[3]] <- male_ori
        data[[rindex]][[4]] <- female_ori

      } else{
        data[[rindex]][[3]] <- matrix(decode(population$breeding[[animals[1]]][[animals[2]]][[index]][[5]]), byrow=FALSE, ncol=4)
        data[[rindex]][[4]] <- matrix(decode(population$breeding[[animals[1]]][[animals[2]]][[index]][[6]]), byrow=FALSE, ncol=4)

      }

      data[[rindex]][[5]] <- colnamed


      rindex <- rindex + 1
    }
  }


  return(data)
}
