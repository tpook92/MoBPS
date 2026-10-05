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

#' Function to compute kinship matrix based on SNPs
#'
#' Function to derive founder haplotypes of selected individuals
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param id Individual IDs to search/collect in the database
#' @param use.id Set TRUE to use animal IDs for column/row-names in the output matrix (default: TRUE)
#' @param array Select a subset of markers to consider in the calculation of kinship
#' @examples
#' data(ex_pop)
#' haplo <- get.inbreeding.snp(ex_pop, gen=2)
#' @return Founder Haplotype data for in gen/database/cohorts selected individuals
#' @export



get.inbreeding.snp <-  function(population=NULL, gen=NULL, database=NULL, cohorts=NULL, id = NULL, use.id=TRUE,
                             array = NULL){


  database <- get.database(population, gen = gen, database = database, cohorts = cohorts, id = id)

  founder_haplo = get.founder.haplotype(population, database = database, array = array)
  ids = names(get.id(population, database = database, use.id = use.id))

  n = ncol(founder_haplo) / 2

  tmp1 = founder_haplo[,(1:n)*2-1]
  tmp2 = founder_haplo[,(1:n)*2]




    inbreeding <- numeric(n)


    for(index in 1:n){

      haplo1_1 <- tmp1[, index]
      haplo2_1 <- tmp2[, index]

      # Haplo1 x Haplo2
      inbreeding[index] = inbreeding[index] + sum(haplo1_1 == haplo2_1)

    }




    inbreeding = inbreeding/nrow(tmp1)

  names(inbreeding) =  ids

  return(inbreeding)
}
