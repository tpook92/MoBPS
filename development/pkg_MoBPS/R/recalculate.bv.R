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

#' Recalculate genomic values
#'
#' Function to recalculate genomic values
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param bv.ignore.traits Vector of traits to ignore in the calculation of the genomic value (default: NULL; Only recommended for high number of traits and experienced users!)
#' @examples
#' data(ex_pop)
#' population <- recalculate.bv(ex_pop, gen=2)
#' @return Population list
#' @export
#'
#'
recalculate.bv <- function(population, gen=NULL, database=NULL, cohorts=NULL, bv.ignore.traits=NULL){

  database <- get.database(population, gen=gen, database = database, cohorts=cohorts)
  store.effect.freq <- FALSE
  import.position.calculation <- NULL
  bit.storing <- FALSE
  nbits <- 30

  if(population$info$miraculix){
    if (requireNamespace("miraculix", quietly = TRUE)) {
      codeOriginsU <- miraculix::codeOrigins
      decodeOriginsU <- miraculix::decodeOrigins
    } else{
      codeOriginsU <- codeOriginsR
      decodeOriginsU <- decodeOriginsR
    }
  } else{
    codeOriginsU <- codeOriginsR
    decodeOriginsU <- decodeOriginsR
  }

  if(length(bv.ignore.traits)>0){
    temp123 <- setdiff(population$info$bv.random.activ , bv.ignore.traits)
  } else{
    temp123 <- population$info$bv.random.activ
  }

  for(index2 in 1:nrow(database)){
    index <- database[index2,1]
    sex <- database[index2,2]
    for(nr.animal in database[index2,3]:database[index2,4]){
      activ_bv <- population$info$bv.random.activ
      if(length(activ_bv)>0){
        temp_out <- calculate.bv(population, index, sex, nr.animal,
                                 activ_bv, import.position.calculation=import.position.calculation,
                                 decodeOriginsU=decodeOriginsU,
                                 store.effect.freq=store.effect.freq,
                                 bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE,
                                 bv.ignore.traits=bv.ignore.traits)
        population$breeding[[index]][[6+sex]][activ_bv,nr.animal] <- temp_out[[1]]
        population$breeding[[index]][[sex]][[nr.animal]][[25]] <- length(bv.ignore.traits)==0

        if(length(temp123)>0){
          population$breeding[[index]][[sex]][[nr.animal]][[26]] <- temp123
        }

        if(store.effect.freq){
          if(length(population$info$store.effect.freq) < index || length(population$info$store.effect.freq[[index]])==0){
            colnames(temp_out[[2]]) <- c("Homo0", "Hetero", "Homo1")
            rownames(temp_out[[2]]) <- population$info$snp.name[population$info$effect.p]
            population$info$store.effect.freq[[index]] <- temp_out[[2]]
          } else{
            population$info$store.effect.freq[[index]] <- population$info$store.effect.freq[[index]] + temp_out[[2]]
          }
        }
      }
    }

  }
 return(population)
}
