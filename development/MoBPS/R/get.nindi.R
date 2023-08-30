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

#' Number of generations
#'
#' Function to calculate the number of generations in the population list
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param count.copy Set to TRUE to double count individuals if multiple copies of an individual are included in gen/database/cohorts
#' @return Numeric value
#' @examples
#' data(ex_pop)
#' get.nindi(ex_pop)
#' @export


get.nindi <- function(population, database=NULL, gen=NULL, cohorts= NULL, extended = FALSE, count.copy = FALSE){

  ids = get.id(population, gen = gen, database = database, cohorts = cohorts)

  if(count.copy){
    nindi = length(ids)
  } else{
    nindi = sum(!duplicated(ids))
  }

  if(extended){

    database = get.database(population, gen = gen, database = database, cohorts = cohorts)

    npheno = get.npheno(population, database = database)
    ngeno = get.genotyped(population, database = database)
    ids = get.id(population, database = database)

    n_pheno = sum(!duplicated(ids[npheno>0]))
    n_geno = sum(!duplicated(ids[ngeno>0]))
    n_both = sum(!duplicated(ids[ngeno>0 & npheno>0]))

    out = c(nindi, n_pheno, n_geno, n_both)
    names(out) = c("total", "phenotyped", "genotyped", "both")
    return(out)
  }

  return(nindi)
}
