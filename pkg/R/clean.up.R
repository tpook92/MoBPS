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

#' Clean-up recombination points
#'
#' Function to remove recombination points + origins with no influence on markers
#' @param population Population list
#' @param gen Generations to clean up (default: "current")
#' @param database Groups of individuals to consider
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @examples
#' data(ex_pop)
#' ex_pop <- clean.up(ex_pop)
#' @return Population-list with deleted irrelevant recombination points
#' @export

clean.up <- function(population, gen="all", database=NULL, cohorts=NULL){
  #remove.no.snp.segments nur auf die Schnelle und super ineffizient!
  generations <- gen
  if(length(gen)==1 && gen=="current"){
    generations <- length(population$breeding)
  }
  if(length(gen)==1 &&  gen=="all"){
    generations <- 1:length(population$breeding)
  }

  database <- get.database(population, generations, database, cohorts)

  for(indexrow in 1:nrow(database)){
    index <- database[indexrow,1]
    index2 <- database[indexrow,2]
    if(length(population$breeding[[index]][[index2]])>0){
      for(index3 in database[indexrow,3]:database[indexrow,4]){
        for(index4 in 1:2){
          removes <- which(diff(population$breeding[[index]][[index2]][[index3]][[index4+4]])==0)+1
          if(length(removes)>0){
            population$breeding[[index]][[index2]][[index3]][[index4+4]] <- population$breeding[[index]][[index2]][[index3]][[index4+4]][-removes]
            population$breeding[[index]][[index2]][[index3]][[index4]] <- population$breeding[[index]][[index2]][[index3]][[index4]][-removes]
          }
        }
      }
    }

  }


  return(population)
}
