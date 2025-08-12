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

#' Set differences between founder pool
#'
#' Function to scale trait for genetic differences based on founder pools
#' @param population Population list
#' @param database Groups of individuals to consider for the export (THIS CAN ONLY BE APPLIED ON FOUNDERS)
#' @param gen Quick-insert for database (vector of all generations to export) (THIS CAN ONLY BE APPLIED ON FOUNDERS, if empty -> default: 1)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export) (THIS CAN ONLY BE APPLIED ON FOUNDERS)
#' @param trait Which trait to set the new mean for
#' @param mean Vector with the target mean for the different pools
#' @param pool Vector with pools considered (default: 1:length(mean))
#' @param reference Target mean is compared again the reference (default: "pool" - average genomic value in the respective pool, alt: "all")
#' @param max.effects Maximum number of locations in the genome that will be assigned an effect for pool-based correction
#' @examples
#' population = creating.diploid(nsnp = 100, nindi = 10, n.additive = 100, founder.pool = 1)
#' population = creating.diploid(population=population, nindi = 10,
#'   founder.pool = 2)
#' population = set.mean.pool(population, mean = c(100,110))
#' @return Class of in gen/database/cohorts selected individuals
#' @export

set.mean.pool = function(population, pool = NULL, mean = NULL, trait = NULL,
                         gen = NULL, database = NULL, cohorts = NULL,
                         reference = "pool",
                         max.effects = Inf){

  if(length(trait)== 0){
    trait = 1
  }

  if(length(pool) == 0){
    pool = 1:length(mean)
  }

  if(length(gen)==0 & length(database)==0 & length(cohorts)==0){
    gen = 1
  }

  database = get.database(population, gen=gen, database = database, cohorts = cohorts)

  founder.pool = get.pool.founder(population, database = database)

  bvs = get.bv(population, database = database)

  if(reference == "pool"){
    mean_per_pool = numeric(length(pool))
    for(index in 1:length(pool)){
      mean_per_pool[index] = mean(bvs[trait, founder.pool==pool[index]])
    }
  } else if(reference == "all"){
    mean_per_pool = rep(mean(bvs[trait, ]), length(pool))

  }


  to_change = mean - mean_per_pool

  chr.nr <- snp.nr <- numeric(sum(population$info$snp))
  nsnp = length(chr.nr)
  pos = 1:nsnp
  start <- 1
  for(index in 1:length(population$info$snp)){
    if(population$info$snp[index]>0){
      chr.nr[start:(start+population$info$snp[index]-1)] <- index
      snp.nr[start:(start+population$info$snp[index]-1)] <- 1:population$info$snp[index]
      start <- start + population$info$snp[index]
    }
  }

  if(max.effects < length(chr.nr)){
    keep1 = sort(sample(length(snp.nr), max.effects))

    snp.nr = snp.nr[keep1]
    chr.nr = chr.nr[keep1]
    nsnp = max.effects
  }
  add_real.bv.add = NULL
  for(index in 1:length(pool)){
    add_real.bv.add = rbind(add_real.bv.add, cbind(snp.nr, chr.nr, 0, to_change[index]/2/ nsnp, to_change[index]/nsnp, pos, pool[index] , TRUE))
  }

  population$info$real.bv.add[[trait]] = rbind(population$info$real.bv.add[[trait]], add_real.bv.add)

  if(population$info$bv.calculated && get.ngen(population)==1){

    for(index5 in 1:nrow(database)){
      temp1 = get.pool.founder(population, database = database[index5,])

      for(indexp in 1:length(pool)){
        population$breeding[[database[index5,1]]][[database[index5,2]+6]][trait,temp1==pool[indexp]] = population$breeding[[database[index5,1]]][[database[index5,2]+6]][trait,temp1==pool[indexp]] + to_change[indexp]
      }
    }
    population$info$pool_effects = TRUE
    population$info$pool_effects_calc = FALSE
  } else{
    population$info$bv.calculated = FALSE
    population$info$bv.calculated.partly = population$info$bv.calculated.partly[population$info$bv.calculated.partly!=trait]
    population$info$pool_effects = TRUE
    population$info$pool_effects_calc = FALSE
    population = breeding.diploid(population)
  }


  return(population)
}

