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
#' @examples
#' data(ex_pop)
#' population <- recalculate.manual(ex_pop, gen = 1)
#' @return Population list
#' @export
#'
#'


recalculate.manual = function(population, gen = NULL, database=NULL, cohorts = NULL, e0 = NULL ,
                              e1 = NULL, e2 = NULL, store.comp.times = TRUE){


  temp123 <- population$info$bv.random.activ

  if(length(temp123) == 0){
    return(population)
  }
  if(store.comp.times){
    tick <- as.numeric(Sys.time())
  }

  if((length(e0)==0 && length(population$info$e0_mat)==0) ||
     (length(e1)==0 && length(population$info$e1_mat)==0) ||
     (length(e2)==0 && length(population$info$e2_mat)==0)){

    effect_matrix0 = effect_matrix1 = effect_matrix2 = matrix(0, ncol = sum(population$info$snp), nrow = population$info$bv.nr)

    for(index in 1:population$info$bv.nr){
      effect_matrix0[index, population$info$real.bv.add[[index]][,6]] = population$info$real.bv.add[[index]][,3]
      effect_matrix1[index, population$info$real.bv.add[[index]][,6]] = population$info$real.bv.add[[index]][,4]
      effect_matrix2[index,population$info$real.bv.add[[index]][,6]] = population$info$real.bv.add[[index]][,5]
    }

    population$info$e0_mat = effect_matrix0
    population$info$e1_mat = effect_matrix1
    population$info$e2_mat = effect_matrix2

    population$info$e0_activ = which(colSums(abs(effect_matrix0))>0)
    population$info$e1_activ = which(colSums(abs(effect_matrix1))>0)
    population$info$e2_activ = which(colSums(abs(effect_matrix2))>0)
  }

  if(length(e0)==0){
    e0  =population$info$e0_mat
  }
  if(length(e1)==0){
    e1  =population$info$e1_mat
  }
  if(length(e2)==0){
    e2  =population$info$e2_mat
  }
  e0_activ = population$info$e0_activ
  e1_activ = population$info$e1_activ
  e2_activ = population$info$e2_activ

  database = get.database(population, gen = gen, database = database, cohorts = cohorts)

  bv_tmp = population$info$bv.ignore.traits.ever

  for(index in 1:nrow(database)){
    geno = get.geno(population, database = database[index,])
    bvs = 0
    if(length(e0_activ)>0){
      g0_tmp = geno[e0_activ,,drop = FALSE]
      #bvs = bvs + (e0)[,e0_activ,drop = FALSE] %*% (g0_tmp == 0L)

      bvs = bvs + (e0)[,e0_activ,drop = FALSE] %*% equal_double(g0_tmp, 0L)
    }
    if(length(e1_activ)>0){
      if(length(e0_activ)>0 && length(e0_activ)==length(e1_activ) && prod(e0_activ==e1_activ)){
        g1_tmp = g0_tmp
      } else{
        g1_tmp = geno[e1_activ,,drop = FALSE]
      }
      #bvs = bvs + (e1)[,e1_activ,drop = FALSE] %*% (g1_tmp==1L)

      bvs = bvs + (e1)[,e1_activ,drop = FALSE] %*% equal_double(g1_tmp,1L)
    }
    if(length(e2_activ)>0){

      if(length(e0_activ)>0 && length(e0_activ)==length(e2_activ) && prod(e0_activ==e2_activ)){
        g2_tmp = g0_tmp
      } else if(length(e1_activ)>0 && length(e1_activ)==length(e2_activ) && prod(e1_activ==e2_activ)){
        g2_tmp = g1_tmp
      } else {
        g2_tmp = geno[e2_activ,,drop = FALSE]
      }
      #bvs = bvs + (e2)[,e2_activ,drop = FALSE] %*% (g2_tmp==2L)

      bvs = bvs + (e2)[,e2_activ,drop = FALSE] %*% equal_double(g2_tmp,2L)
    }
    bvs = bvs + population$info$base.bv

    population$breeding[[database[index,1]]][[database[index,2]+6]][,database[index,3]:database[index,4]] = bvs

    if(bv_tmp){
      for(index2 in database[index,3]:database[index,4]){
        population$breeding[[database[index,1]]][[database[index,2]]][[index2]][[25]] <- TRUE
        population$breeding[[database[index,1]]][[database[index,2]]][[index2]][[26]] <- temp123
      }
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


