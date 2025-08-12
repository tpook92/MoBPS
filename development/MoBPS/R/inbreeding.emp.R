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

#' Empirical kinship
#'
#' Function to compute empirical kinship for a set of individuals)
#' @param population Population list
#' @param animals List of animals to compute kinship for
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param use.id Set TRUE to use animal IDs for column/row-names in the output matrix (default: TRUE)
#' @examples
#' data(ex_pop)
#' inbreeding <- inbreeding.emp(population=ex_pop, database=cbind(3,1,1,25))
#' @return Empirical kinship matrix (IBD-based since Founders)
#' @export

inbreeding.emp <- function(population=NULL, animals=NULL, gen=NULL, database=NULL, cohorts=NULL, use.id=TRUE){


  if(length(animals)==0 || length(gen)>0 || length(database)>0 || length(cohorts)>0){
    database <- get.database(population, gen, database, cohorts)
    animals <- list()
    for(index in 1:nrow(database)){
      if(diff(database[index,3:4])>=0){
        for(index2 in database[index,3]:database[index,4]){
          animals[[length(animals)+1]] <- population$breeding[[database[index,1]]][[database[[index,2]]]][[index2]]
        }
      }

    }
  }
  n <- length(animals)
  inbreeding <- numeric(n)

  col_names <- get.id(population, database=database)
  if(use.id==FALSE){
    col_names <- names(col_names)
  }


  chrom.length <- max(animals[[1]][[1]])
  for( i in 1:n){
    j = i
      chr <- list()
      chr[[1]] <- animals[[i]][[1]][-1]
      chr[[2]] <- animals[[i]][[2]][-1]
      chr[[3]] <- animals[[j]][[1]][-1]
      chr[[4]] <- animals[[j]][[2]][-1]
      origin <- list()
      origin[[1]] <- animals[[i]][[5]]
      origin[[2]] <- animals[[i]][[6]]
      origin[[3]] <- animals[[j]][[5]]
      origin[[4]] <- animals[[j]][[6]]

      activ <- c(1,1,1,1)
      prev <- 0
      activ.recom <- c(chr[[1]][activ[1]], chr[[2]][activ[2]], chr[[3]][activ[3]], chr[[4]][activ[4]])
      activ.ursprung <- c(origin[[1]][activ[1]],origin[[2]][activ[2]],origin[[3]][activ[3]],origin[[4]][activ[4]])



      for(steps in 1:(length(c(chr[[1]], chr[[2]], chr[[3]], chr[[4]]))-3)){
        activ.min <- which.min(activ.recom)[1]
        activ.posi <- chr[[activ.min]][activ[activ.min]]

        ibd.factor <- (sum(activ.ursprung[1] == activ.ursprung[3:4]) + sum(activ.ursprung[2] == activ.ursprung[3:4]))/4
        inbreeding[i] <- inbreeding[i] + ibd.factor * (activ.posi - prev) / chrom.length
        prev <- activ.posi

        activ[activ.min] <- min(activ[activ.min] +1, length(chr[[activ.min]]))
        activ.recom[activ.min] <- chr[[activ.min]][activ[activ.min]]
        activ.ursprung[activ.min] <-  origin[[activ.min]][activ[activ.min]]

      }
    }


  names(inbreeding) <- col_names

  inbreeding = (inbreeding - 0.5) * 2

  return(inbreeding)
}
