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

#' Mutation intro
#'
#' Function to change the base-pair in a specific loci
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param qtl.posi Marker number to mutate
#' @param target.variant target variant to obtain (( if haplotype already is correct do not introduce a mutation ))
#' @param haplo.set Select chromosome set (default: 1 , alt: 2, 1:2 (to edit both))
#' @examples
#' data(ex_pop)
#' ex_pop <- mutation.intro(ex_pop, database = cbind(1,1,1), qtl.posi=100)
#' @return Population-List with mutated marker for the selected individual
#' @export

mutation.intro <- function(population, gen = NULL, database = NULL, cohorts = NULL, qtl.posi, target.variant = NULL, haplo.set=1) {

  database <- get.database(population, gen, database, cohorts)

  for(row in 1:nrow(database)){
    for(index in database[row,3]:database[row,4]){

      gen = database[row,1]
      sex = database[row,2]
      individual.nr = index

      if(length(target.variant)>0){
        current_haplo = get.haplo(population, database = cbind(database[row,1], database[row,2], index,index))[qtl.posi,]
      } else{
        current_haplo = c(-1,-1) #
      }
      for(activ.set in haplo.set){
        if(length(target.variant)==0 || current_haplo[activ.set] != target.variant){

          if(sum(population$breeding[[gen]][[sex]][[individual.nr]][[2+ activ.set]]==qtl.posi)==0){
            population$breeding[[gen]][[sex]][[individual.nr]][[2+ activ.set]] <- sort(c(qtl.posi,population$breeding[[gen]][[sex]][[individual.nr]][[2+ activ.set]]))
          } else{
            population$breeding[[gen]][[sex]][[individual.nr]][[2+ activ.set]] <- unique(c(qtl.posi,population$breeding[[gen]][[sex]][[individual.nr]][[2+ activ.set]]))[-1]
          }
        }

      }


    }
  }

  return(population)
}
