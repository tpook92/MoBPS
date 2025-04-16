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

#' Function to exclude individuals from a database
#'
#' Function to exclude individuals from a database
#' @param population Population list
#' @param old.name Quick-insert for database (vector of names of cohorts to export)
#' @param new.name Groups of individuals to remove from the database (same IDs!)
#' @return population list
#' @examples
#' data(ex_pop)
#' population <- rename.cohort(ex_pop, old.name="Cohort_4", new.name = "NewName")
#' @export

rename.cohort <- function(population, old.name, new.name, verbose = TRUE){

  to_rename = which(population$info$cohorts[,1] == old.name)


  if(length(to_rename) == 0){

    old.name1 = paste0(old.name, "_M")

    to_rename = which(population$info$cohorts[,1] == old.name1)
    if(length(to_rename)>0){
      population$info$cohorts[to_rename,1] = paste0(new.name, "_M")
      rownames(population$info$cohorts)[to_rename] = paste0(new.name, "_M")
      if(verbose){cat(paste0("Cohort ",old.name1, " as been renamed to ", paste0(new.name, "_M"), "\n"))}
    }



    old.name2 = paste0(old.name, "_F")

    to_rename = which(population$info$cohorts[,1] == old.name2)
    if(length(to_rename)>0){
      population$info$cohorts[to_rename,1] = paste0(new.name, "_F")
      rownames(population$info$cohorts)[to_rename] = paste0(new.name, "_F")
      if(verbose){cat(paste0("Cohort ",old.name2, " as been renamed to ", paste0(new.name, "_F"), "\n"))}

    }
  } else{
    population$info$cohorts[to_rename,1] = new.name
    rownames(population$info$cohorts)[to_rename] = new.name
    if(verbose){cat(paste0("Cohort ",old.name, " as been renamed to ", new.name, "\n"))}

  }

  return(population)

}
