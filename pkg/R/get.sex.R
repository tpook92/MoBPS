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

#' Extraction of individual sex
#'
#' Function to extract the sex of selected individuals
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names
#' @param male.female.coding Set to TRUE to display male/female instead of 1/2
#' @export
#'


get.sex <- function(population, database=NULL, gen=NULL, cohorts=NULL, use.id=F, male.female.coding = F){
  database <- get.database(population, gen, database, cohorts)
  sex <- get.full.database(population = population, database = database)[,2]
  ids <- get.id(population = population, database = database, use.id = use.id)
  names(sex)<- ids

  if(male.female.coding){
    sex[sex == 1] <- "male"
    sex[sex == 2] <- "female"
  }

  return(sex)
}
