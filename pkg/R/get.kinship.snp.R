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
#' @param mode Set to 2 to not use C++ based implementation (Why would you ;) )
#' @param array Select a subset of markers to consider in the calculation of kinship
#' @examples
#' data(ex_pop)
#' haplo <- get.kinship.snp(ex_pop, gen=2)
#' @return Founder Haplotype data for in gen/database/cohorts selected individuals
#' @export



get.kinship.snp <-  function(population=NULL, gen=NULL, database=NULL, cohorts=NULL, id = NULL, use.id=TRUE,
                             mode = 1, array = NULL){


  database <- get.database(population, gen = gen, database = database, cohorts = cohorts, id = id)

  founder_haplo = get.founder.haplotype(population, database = database, array = array)
  ids = names(get.id(population, database = database, use.id = use.id))

  n = ncol(founder_haplo) / 2




  tmp1 = founder_haplo[,(1:n)*2-1]
  tmp2 = founder_haplo[,(1:n)*2]


  if(mode == 1){

    kinship = kinship_cpp(tmp1, tmp2)

  } else{

    kinship <- matrix(0, nrow=n, ncol=n)


    for(index in 1:n){

      haplo1_1 <- tmp1[, index]
      haplo2_1 <- tmp2[, index]

      # Haplo1 x Haplo1
      kinship[,index] = kinship[,index] + colSums(tmp1 == haplo1_1)
      # Haplo1 x Haplo2
      kinship[,index] = kinship[,index] + colSums(tmp1 == haplo2_1)

      # Haplo2 x Haplo1
      kinship[,index] = kinship[,index] + colSums(tmp2 == haplo1_1)
      # Haplo2 x Haplo2
      kinship[,index] = kinship[,index] + colSums(tmp2== haplo2_1)

    }
  }



  kinship = kinship/4/nrow(tmp1)

  colnames(kinship) = rownames(kinship) = ids

  return(kinship)
}


if(FALSE){

  library(Rcpp)
  sourceCpp("kinship.cpp")


  # just some testing if its faster
  population = creating.diploid(nsnp = 10000, nindi = 100, chr.nr = 5, nqtl = 100, n.additive = 100)
  for(index in 1:10){
    population = breeding.diploid(population, breeding.size = c(1000,0))
  }
  population = clean.up(population, gen = 11)
  system.time({A1 = get.kinship.emp(population, database = cbind(11,1,1,200))})
  system.time({A2 = get.kinship.snp(population, database = cbind(11,1,1,200), mode = 1)})
  system.time({A3 = get.kinship.snp(population, database = cbind(11,1,1,200), mode = 2)})
  system.time({A4 = get.kinship.snp(population, database = cbind(11,1,1,200), mode = 1, array = "onlyQTL")})
  system.time({A5 = get.kinship.snp(population, database = cbind(11,1,1,200), mode = 2, array = "onlyQTL")})

  max(diff(A2-A1)) # inaccuracy recombination to snps
  max(diff(A3-A2)) # check
  max(diff(A5-A3)) # inaccuracy density reduction

}
