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

#' Derive genotypes of selected individuals
#'
#' Function to derive genotypes of selected individuals
#' @param population Population list
#' @param geno Manually provided genotype matrix (SNP x Indi)
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param id Individual IDs to search/collect in the database
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names (default: TRUE)
#' @param array Use only markers available on the array
#' @param p.database Groups of individuals to consider for the export
#' @param p.gen Quick-insert for database (vector of all generations to export)
#' @param p.cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param p.id Individual IDs to search/collect in the database
#' @examples
#' data(ex_pop)
#' G <- get.grm(ex_pop, gen=2)
#' @return Genotype data for in gen/database/cohorts selected individuals
#' @export

get.grm = function(population, geno = NULL, database=NULL, gen=NULL, cohorts=NULL, id = NULL, use.id=TRUE, array = NULL,
                   p.database = NULL, p.gen = NULL, p.cohorts = NULL, p.id = NULL){


  if(length(geno)==0){
    database = get.database(population, database = database, gen = gen, cohorts = cohorts, id = id)
    geno = get.geno(population, database = database, use.id = use.id, array = array)
  }

  if(length(p.database)>0 || length(p.gen) > 0 || length(p.cohorts)>0 || length(p.id)>0){
    p.database = get.database(population, database = p.database, gen = p.gen, cohorts = p.cohorts, id = p.id)

    p.geno = get.geno(population, database = p.database, use.id = use.id, array = array)

    p = rowMeans(p.geno)/2
  } else{
    p = rowMeans(geno)/2
  }

  if(population$info$miraculix){
    Z.code = miraculix::genomicmatrix(geno)
    G <- miraculix::relationshipMatrix(Z.code, centered=FALSE, normalized=FALSE)
    G <- scaling.relationship(G, Z.code, p)
    colnames(G) = rownames(G) = colnames(geno)
  } else{
    Ztm <- geno - p * 2
    G <- crossprod(Ztm)/ (2 * sum(p*(1-p)))

  }

  return(G)

}
