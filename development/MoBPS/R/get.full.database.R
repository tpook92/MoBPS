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

#' Generation of a database with one row per individual
#'
#' Function to generate a database with one row per individual
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param id Individual IDs to search/collect in the database
#' @export
#'
#'
get.full.database <- function(population, database=NULL, gen=NULL,
                              cohorts=NULL, id=NULL){
  database <- get.database(population, gen = gen, database = database,
                           cohorts = cohorts, id=id)

  out.db <- c()
  for(i in 1:nrow(database)){
    if(database[i, 3] == database[i, 4]){
      out.db <- rbind(out.db, database[i,])
    } else{
      animal_numbers <- database[i, 3]:database[i, 4]
      new_chunk <- cbind(matrix(database[i,1:2], ncol = 2,
                                nrow = length(animal_numbers),
                                byrow = T), animal_numbers, animal_numbers)
      out.db <- rbind(out.db, new_chunk)
    }
  }

  return(out.db)
}
