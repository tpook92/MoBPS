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

#' Merging of cohorts
#'
#' Function to merge cohorts in a population list
#' @param population Population list
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param name.cohort Name of the newly added cohort
#' @examples
#' data(ex_pop)
#' population <- breeding.diploid(ex_pop, add.gen = 7, breeding.size = 50)
#' population <- breeding.diploid(population, add.gen = 7, breeding.size = 50,
#'  selection.m.database = cbind(6,1), selection.f.database = cbind(6,2))
#' population <- merging.cohorts(population, cohorts = c("Cohort_7_M", "Cohort_8_M"))
#' @export


merging.cohorts = function(population, cohorts, name.cohort = NULL){

  cohort_list = get.cohorts(population, extended = TRUE)

  to_merge = length(cohorts)
  for(index in 1:length(cohorts)){
    to_merge[index] = which(cohort_list[,1]==cohorts[index])
  }

  merging = cohort_list[to_merge,]
  suppressWarnings(  storage.mode(merging) <- "numeric")


  if(length(name.cohort)==0){
    name.cohort <- paste0("Cohort_", population$info$cohort.index)
    population$info$cohort.index <- population$info$cohort.index + 1
  }


  new_cohort = c(name.cohort, mean(merging[,2]), sum(merging[,3]), sum(merging[,4]), mean(merging[,5]),
                 min(merging[,6]), min(merging[,7]), mean(merging[,8]), mean(merging[,9]),
                 min(as.numeric(merging[,10])), max(as.numeric(merging[,11])))

  # some checks to make sure what you are doing is ok

  if(sum(new_cohort[3:4]>0)>1){
    stop("Illegal merge! Only cohorts of the same sex can be merged")
  }

  if(sum(new_cohort[2] == merging[,2])< nrow(merging)){
    stop("Illegal merge! Only cohorts of the same generation can be merged")
  }

  if(sum(new_cohort[5] == merging[,5])< nrow(merging)){
    stop("Illegal merge! Only cohorts of the same class can be merged")
  }

  if(new_cohort[4]>0){
    if((min(merging[,7])+as.numeric(new_cohort[4])) != (max(merging[,7]) + merging[which.max(merging[,7]),4])){
      stop("Illegal merge! Only adjacent cohorts can be merged")
    }
  }

  if(new_cohort[3]>0){
    if((min(merging[,6])+as.numeric(new_cohort[3])) != (max(merging[,6]) + merging[which.max(merging[,6]),3])){
      stop("Illegal merge! Only adjacent cohorts can be merged")
    }
  }

  population$info$cohorts = rbind(cohort_list[-to_merge,], new_cohort)

  return(population)
}
