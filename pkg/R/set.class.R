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
#' data(ex_pop)
#' population <- set.class(ex_pop, database=cbind(1,1), new.class = 2)
#' @return Population-List with newly entered class values
#' @export

set.class <- function(population, database=NULL, gen=NULL, cohorts=NULL, new.class=0){

  database <- get.database(population, gen, database, cohorts)
  for(index in 1:nrow(database)){
    population$breeding[[database[index,1]]][[database[index,2]+4]][database[index,3]:database[index,4]] <- new.class
  }
  return(population)

}
