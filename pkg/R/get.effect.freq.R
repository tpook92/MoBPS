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

#' Compute marker frequency in QTL-markers
#'
#' Function to compute marker frequency in QTL-markers
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param sort Set to FALSE to not sort markers according to position on the genome
#' @examples
#' data(ex_pop)
#' get.effect.freq(ex_pop, gen=1)
#' @return Matrix with allele frequencies in the QTLs
#' @export

get.effect.freq <- function(population, database=NULL, gen=NULL, cohorts=NULL, sort=FALSE){

  order <- population$info$effect.p
  if(sort){
    order <- sort(order)
  }
  genos <- get.geno(population, gen = gen, database = database, cohorts = cohorts)[order,]


  effects <- cbind(rowSums(genos==0), rowSums(genos==1), rowSums(genos==2))
  colnames(effects) <- c("Homo0", "Hetero", "Homo1")

  return(effects)

}
