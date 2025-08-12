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

#' Expected inbreeding
#'
#' Function to derive pedigree based inbreeding
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param depth.pedigree Depth of the pedigree in generations
#' @param start.kinship Relationship matrix of the individuals in the first considered generation
#' @param elements Vector of individuals from the database to include in pedigree matrix
#' @param storage.save Lower numbers will lead to less memory but slightly higher computing time (default: 1.5, min: 1)
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' data(ex_pop)
#' inbreeding <- inbreeding.exp(population=ex_pop, gen=5)
#' @return Pedigree-based inbreeding in gen/database/cohort selected individuals
#' @export
#'

inbreeding.exp <- function(population, gen=NULL, database=NULL, cohorts=NULL, depth.pedigree=7,
                           start.kinship=NULL,
                           elements = NULL,
                           storage.save=1.5,
                           verbose=TRUE){

  A <- kinship.exp(population = population, gen=gen, database=database,
                   cohorts=cohorts, depth.pedigree=depth.pedigree,
                   start.kinship=start.kinship,
                   elements=elements,
                   storage.save=storage.save,
                   verbose=verbose,
                   mult=2)

  inbreeding = (diag(A)-1)
  return(inbreeding)
}
