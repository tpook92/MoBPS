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

#' Recalculate genomic values
#'
#' Function to recalculate genomic values
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param e0 Effect matrix for 0 genotype (default: Will be automatically extracted)
#' @param e1 Effect matrix for 1 genotype (default: Will be automatically extracted)
#' @param e2 Effect matrix for 2 genotype (default: Will be automatically extracted)
#' @param store.comp.times If TRUE store computation times in $info$comp.times.general (default: TRUE)
#' @return Population list
#' @export
#'
#'


recalculate.manual = function(population, gen = NULL, database=NULL, cohorts = NULL, e0 = NULL ,
                              e1 = NULL, e2 = NULL, store.comp.times = TRUE){


  temp123 <- population$info$bv.random.activ

  if(store.comp.times){
    tick <- as.numeric(Sys.time())
  }

  if((length(e0)==0 && length(population$info$e0)==0) ||
     (length(e1)==0 && length(population$info$e1)==0) ||
     (length(e2)==0 && length(population$info$e2)==0)){

    effect_matrix0 = effect_matrix1 = effect_matrix2 = matrix(0, ncol = sum(population$info$snp), nrow = population$info$bv.nr)

    for(index in 1:population$info$bv.nr){
      effect_matrix0[index, population$info$real.bv.add[[index]][,6]] = population$info$real.bv.add[[index]][,3]
      effect_matrix1[index, population$info$real.bv.add[[index]][,6]] = population$info$real.bv.add[[index]][,4]
      effect_matrix2[index,population$info$real.bv.add[[index]][,6]] = population$info$real.bv.add[[index]][,5]
    }

    population$info$e0 = effect_matrix0
    population$info$e1 = effect_matrix1
    population$info$e2 = effect_matrix2
  }

  if(length(e0)==0){
    e0  =population$info$e0
  }
  if(length(e1)==0){
    e1  =population$info$e1
  }
  if(length(e2)==0){
    e2  =population$info$e2
  }

  database = get.database(population, gen = gen, database = database, cohorts = cohorts)

  for(index in 1:nrow(database)){
    geno = get.geno(population, database = database[index,])
    bvs = (e0) %*% (geno==0) + (e1) %*% (geno==1) + (e2) %*% (geno==2) + population$info$base.bv
    population$breeding[[database[index,1]]][[database[index,2]+6]][,database[index,3]:database[index,4]] = bvs

    for(index2 in database[index,3]:database[index,4]){
      population$breeding[[database[index,1]]][[database[index,2]]][[index2]][[25]] <- TRUE
      population$breeding[[database[index,1]]][[database[index,2]]][[index2]][[26]] <- temp123

    }
  }

  if(store.comp.times){
    tock <- as.numeric(Sys.time())


    comp.times <- c(0, tock - tick, 0,0,0,0, tock - tick)
    comp.times[comp.times<0] <- 0
    comp.times[comp.times>10e6] <- 0

    population$info$comp.times.general<- round(rbind(population$info$comp.times.general, comp.times, deparse.level = 0), digits=4)
    if(nrow(population$info$comp.times.general)==1){
      colnames(population$info$comp.times.general) <- c("preparation", "new real BV", "phenotypes", "BVE","selection","generate new individuals","total")
    }
  }


  return(population)
}


