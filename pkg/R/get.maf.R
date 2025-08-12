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

#' Calculate minor allele frequencies
#'
#' Function to calculate minor allele frequencies
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @examples
#' data(ex_pop)
#' get.maf(ex_pop, gen = 1)
#' @return Allele frequency of the minor allele
#' @export

get.maf <- function(population, database=NULL, gen=NULL, cohorts=NULL){

  database <- get.database(population, gen, database, cohorts)

  p_i = get.allele.freq(population, database = database)

  p_i[p_i>0.5] = 1 - p_i[p_i>0.5]

  return(p_i)

}
