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

#' Export underlying phenotypes
#'
#' Function to export underlying phenotypes
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param use.all.copy Set to TRUE to extract phenotyping
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names (default: FALSE)
#' @examples
#' data(ex_pop)
#' get.pheno(ex_pop, gen=2)
#' @return Phenotypes for in gen/database/cohorts selected individuals
#' @export

get.pheno.single <- function(population, database=NULL, gen=NULL, cohorts=NULL, use.all.copy = FALSE, use.id=FALSE){

  database <- get.database(population, gen, database, cohorts)

  n.animals <- sum(database[,4] - database[,3] +1)
  data <- list()
  before <- 1
  names <- numeric(n.animals)

  if(use.all.copy){

    for(row in 1:nrow(database)){
      animals <- database[row,]

      for(index in database[row,3]:database[row,4]){
        copies <- population$breeding[[animals[1]]][[animals[2]]][[index]][[21]]
        n_obs <- numeric(population$info$bv.nr)

        phenos <- matrix(NA, nrow=population$info$bv.nr, ncol = 0)

        for(rows in 1:nrow(copies)){

          if(sum(population$breeding[[copies[rows,1]]][[copies[rows,2]]][[copies[rows,3]]][[15]]> n_obs)>=1){
            switches <- population$breeding[[copies[rows,1]]][[copies[rows,2]]][[copies[rows,3]]][[15]]> n_obs

            if(sum(switches)>0){
              pheno.indi <- population$breeding[[copies[rows,1]]][[copies[rows,2]]][[copies[rows,3]]][[27]]
            }
            n_obs[switches] <- population$breeding[[copies[rows,1]]][[copies[rows,2]]][[copies[rows,3]]][[15]][switches]

            if(max(n_obs)>ncol(phenos)){
              phenos <- cbind(phenos, matrix(NA, nrow=population$info$bv.nr, ncol = max(n_obs) - ncol(phenos)))
            }

            for(bven in which(switches)){
              if(n_obs[bven]>0){
                data[[before]][bven,1:n_obs[bven]] <- pheno.indi[[bven]]
              }

            }


          }

        }
        names[before] <- paste(if(animals[2]==1) "M" else "F", index,"_", animals[1], sep="")

        before <- before +1

      }
    }
  } else{
    for(row in 1:nrow(database)){

      animals <- database[row,]
      for(nr.animal in animals[3]:animals[4]){

        n.obs <- (population$breeding[[animals[1]]][[animals[2]]][[nr.animal]][[15]])
        pheno.indi <- (population$breeding[[animals[1]]][[animals[2]]][[nr.animal]][[27]])
        phenos <- matrix(NA, nrow=population$info$bv.nr, ncol = max(n.obs))

        for(bven in 1:population$info$bv.nr){
          if(n.obs[bven]>0){
            phenos[bven, 1:n.obs[bven]] <- pheno.indi[[bven]]
          }

        }

        data[[before]] <- phenos
        names[before] <- paste(if(animals[2]==1) "M" else "F", nr.animal,"_", animals[1], sep="")
        before <- before + 1
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
