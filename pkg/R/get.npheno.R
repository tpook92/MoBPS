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

#' Export underlying number of observations per phenotype
#'
#' Function to export the number of observation of each underlying phenotype
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param use.all.copy Set to TRUE to extract phenotyping
#' @examples
#' data(ex_pop)
#' get.pheno(ex_pop, gen=2)
#' @return Phenotypes for in gen/database/cohorts selected individuals
#' @export

get.npheno <- function(population, database=NULL, gen=NULL, cohorts=NULL, use.all.copy = FALSE){

  database <- get.database(population, gen, database, cohorts)

  n.animals <- sum(database[,4] - database[,3] +1)
  data <- matrix(NA, ncol=n.animals, nrow=population$info$bv.nr)
  before <- 0
  names <- numeric(n.animals)

  if(use.all.copy){

    for(row in 1:nrow(database)){
      animals <- database[row,]
      nanimals <- database[row,4] - database[row,3] +1
      names[(before+1):(before+nanimals)] <- paste(if(animals[2]==1) "M" else "F", animals[3]:animals[4],"_", animals[1], sep="")

      for(index in database[row,3]:database[row,4]){
        copies <- population$breeding[[animals[1]]][[animals[2]]][[index]][[21]]
        n_obs <- numeric(population$info$bv.nr)

        for(rows in 1:nrow(copies)){

          if(sum(population$breeding[[copies[rows,1]]][[copies[rows,2]]][[copies[rows,3]]][[15]]> n_obs)>=1){
            switches <- population$breeding[[copies[rows,1]]][[copies[rows,2]]][[copies[rows,3]]][[15]]> n_obs
            data[switches, (before+1)] <- population$breeding[[copies[rows,1]]][[copies[rows,2]+8]][switches,copies[rows,3]]
            n_obs[switches] <- population$breeding[[copies[rows,1]]][[copies[rows,2]]][[copies[rows,3]]][[15]][switches]
          }

        }

        before <- before +1

      }
    }
  } else{
    for(row in 1:nrow(database)){
      animals <- database[row,]
      nanimals <- database[row,4] - database[row,3] +1
      if(nanimals>0){

        names[(before+1):(before+nanimals)] <- paste(if(animals[2]==1) "M" else "F", animals[3]:animals[4],"_", animals[1], sep="")
        for(index in animals[3]:animals[4]){
          before <- before + 1
          data[,before] <- population$breeding[[animals[1]]][[animals[2]]][[index]][[15]]
        }
      }

    }
  }

  row_names <- paste("Trait", 1:population$info$bv.nr)
  colnames(data) <- names
  rownames(data) <- row_names
  return(data)
}
