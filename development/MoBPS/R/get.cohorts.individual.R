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

#' Derive ID on an individual
#'
#' Function to derive the internal ID given to each individual
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names
#' @param keep.order To not change order of individuals when ids are provided (default: FALSE)
#' @examples
#' data(ex_pop)
#' get.cohorts.individual(ex_pop, gen=2)
#' @return Individual ID for in gen/database/cohorts selected individuals
#' @export

get.cohorts.individual <- function(population, database=NULL, gen=NULL, cohorts=NULL, use.id=FALSE, keep.order = FALSE){

  database <- get.database(population, gen, database, cohorts, keep.order = keep.order, per.individual = TRUE)

  cohorts_list = get.cohorts(population, extended = TRUE)[,c(2,3,4,6,7)]
  storage.mode(cohorts_list) = "numeric"

  cohorts_names = rownames(cohorts_list)

  cohort_included = names = numeric(nrow(database))
  for(index in 1:nrow(database)){

    gen_check = cohorts_list[,1] == database[index,1]
    sex_check = database[index,2] == (1 + (cohorts_list[,3]>0))
    nr_check = database[index,3] >= cohorts_list[,3+database[index,2]] & (database[index,3] < (cohorts_list[,3+database[index,2]] + cohorts_list[,1+database[index,2]]))

    names[index] = paste0(if(database[index,2]==1) "M" else "F", database[index,3],"_", database[index,1])

    cohort_included[index] = cohorts_names[gen_check & sex_check & nr_check]
  }

  if(use.id){
    names(cohort_included) <- get.id(population, database = database)
  } else{
    names(cohort_included) <- names
  }

  return(cohort_included)
}
