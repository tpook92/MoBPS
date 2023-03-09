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

#' Calculate breeding values
#'
#' Internal function to calculate the breeding value of a given individual
#' @param population Population list
#' @param breeding.size Number of individuals to generate (default: 100)
#' @param selection.rate Proportion of individuals to select in each breeding cycle (default: 0.5)
#' @param pool.gen Generations of individuals to consider as founder pool to start from (default: NULL)
#' @param pool.cohorts Cohorts of individuals to consider as founder pool to start from (default: NULL)
#' @param pool.database Groups of individuals to consider as founder pool to start from (default: NULL)
#' @param target.gen Generations of individuals to consider to calculate target genomic value to get to (default: NULL)
#' @param target.cohorts Cohorts of individuals to consider to calculate target genomic value to get to (default: NULL)
#' @param target.database Groups of individuals to consider to calculate target genomic value to get to (default: NULL)
#' @param reduction.multiplier Traits that already exceed the target are bred against. Weighting is scaled by this factor.
#' @param name.cohort Name of the newly added cohort
#' @param sex.quota Share of newly added female individuals (default: 0.5)
#' @param add.gen Generation you want to add the new individuals to (default: 1)
#' @param target.direction Manual select with traits are supposed to increase / decrease (1 target high, -1 target low)
#' @param verbose Set to FALSE to not display any prints
#' @return population list with newly added individuals
#' @export

add.diversity = function(population,
                         breeding.size = 100,
                         selection.rate = 0.5,
                         pool.gen = NULL, pool.database = NULL, pool.cohorts = NULL,
                         target.gen = NULL, target.database = NULL, target.cohorts = NULL,
                         reduction.multiplier = 5,
                         name.cohort = NULL,
                         sex.quota = NULL,
                         add.gen = 1,
                         target.direction = NULL,
                         verbose = TRUE){


  pool.database = get.database(population, gen = pool.gen, database = pool.database, cohorts=  pool.cohorts)
  target.database = get.database(population, gen = target.gen, database = target.database, cohorts = target.cohorts)

  target_bv <- rowMeans(get.bv(population, database = target.database))
  temp1 <- population$info$real.bv.add
  temp2 <- population$info$real.bv.mult
  temp3 <- population$info$real.bv.dice
  temp1[[length(temp1)]] <- NULL
  temp2[[length(temp2)]] <- NULL
  temp3[[length(temp3)]] <- NULL

  dataset = get.haplo(population, database = pool.database)

  map = get.map(population)

  pop1 = creating.diploid(dataset = dataset, map = map, real.bv.add = temp1, real.bv.mult = temp2, real.bv.dice = temp3,
                          verbose = FALSE)

  new_bvs <- rowMeans(get.bv(pop1, gen = 1))

  if(length(target.direction)==0){
    target.direction = 2*((target_bv > new_bvs)-0.5)
  }

  j <- 1

  while(sum(new_bvs*target.direction)<sum(target_bv * target.direction)){

    cur <- length(pop1$breeding)

    # some slow selection to get a population with similar genomic value but different genotypes

    weights = target_bv - new_bvs
    weights[weights<0 & target.direction==1] =   weights[weights<0 & target.direction==1] * reduction.multiplier
    weights[weights>0 & target.direction==(-1)] =   weights[weights>0 & target.direction==(-1)] * reduction.multiplier
    pop1 <- breeding.diploid(pop1,
                             breeding.size = c(sum(breeding.size),0),
                             selection.size = c(sum(breeding.size)* selection.rate,0),
                             selection.criteria = "bv",
                             multiple.bve.weights.m = weights,
                             verbose=FALSE)

    new_bvs <- rowMeans(get.bv(pop1, gen=cur+1))

    j <- j+1
  }

  haplo_new <- get.haplo(pop1, gen= max(1,(j-1)))
  new_bvs <- rowMeans(get.bv(pop1, gen= max(1,(j-1))))

  if(verbose) cat(paste0("Required ", j , " generations to generate introduced material.\n"))
  if(verbose) cat(paste0("Avg. genomic value for traits:\n"))
  if(verbose) print(new_bvs)
  if(verbose) cat(paste0("Compared to target reference:\n"))
  if(verbose) print(target_bv)

  if(length(breeding.size)==2 && length(sex.quota)==0){
    sex.quota = breeding.size[2] / sum(breeding.size)
  } else if(length(sex.quota)==0){
    sex.quota = 0.5
  }

  population <- creating.diploid(population = population,
                                 dataset = haplo_new,
                                 sex.quota = sex.quota,
                                 name.cohort = name.cohort,
                                 generation = add.gen,
                                 verbose = verbose)
}
