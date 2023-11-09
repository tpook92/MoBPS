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

#' Export estimated breeding values
#'
#' Function to export estimated breeding values
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param new.class Class to change to (either single character or vector for each individual when just a single group is selected)
#' @examples
#' data(ex_pop)
#' population <- set.time.point(ex_pop, database=cbind(1,1), time.point = 2)
#' @return Population-List with newly entered class values
#' @export

set.age.point <- function(population, database=NULL, gen=NULL, cohorts=NULL, time.point=0){

  database <- get.database(population, gen, database, cohorts)
  if(length(database)>0){
    if(length(time.point)==1){
      for(index in 1:nrow(database)){
        population$breeding[[database[index,1]]][[database[index,2]+22]][database[index,3]:database[index,4]] <- time.point
      }
    } else{
      next1 = 1
      for(index in 1:nrow(database)){
        n_indi = database[index,4] - database[index,3] + 1
        population$breeding[[database[index,1]]][[database[index,2]+22]][database[index,3]:database[index,4]] <- time.point[next1:(next1 + n_indi - 1)]
        next1 = next1 + n_indi
      }

    }

  }

  return(population)

}

