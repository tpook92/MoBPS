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

#' QTL extraction
#'
#' Function to the position of QTLs (for snp/chr use get.qtl.effects()
#' @param population Population list
#' @examples
#' data(ex_pop)
#' positions <- get.qtl(ex_pop)
#' @return Vector of SNP positions
#' @export

get.qtl <- function(population){

  return(population$info$effect.p)

}


#' QTL effect extraction
#'
#' Function to extract QTL effect sizes
#' @param population Population list
#' @examples
#' data(ex_pop)
#' effects <- get.qtl.effects(ex_pop)
#' @return List with [[1]] single SNP QTLs [[2]] epistatic SNP QTLs [[3]] dice QTL
#' @export

get.qtl.effects <- function(population){

  adds <- population$info$real.bv.add
  mults <- population$info$real.bv.mult
  dice <- population$info$real.bv.dice

  adds <- adds[-length(adds)]
  mults <- mults[-length(mults)]

  for(index in 1:length(adds)){
    if(length(adds[[index]])>0){
      colnames(adds[[index]]) <- c("SNP", "Chromo", "Effect AA", "Effect AB", "Effect BB", "Overall position", "Pool")
    }
  }
  for(index in 1:length(mults)){
    if(length(mults[[index]])>0){
      colnames(mults[[index]]) <- c("SNP 1", "Chromo 1", "SNP 2", "Chromo 2", "Effect 00", "Effect 01", "Effect 02", "Effect 10", "Effect 11", "Effect 12",
                           "Effect 20", "Effect 21", "Effect 22")
    }
  }
  dice <- dice[-length(dice)]
  effects <- list(adds, mults, dice)

  return(effects)

}

#' QTL effect variance extraction
#'
#' Function to extract QTL effect variance for single SNP QTLs in a given gen/database/cohort
#' @param population Population list
#' @param database Groups of individuals to consider
#' @param gen Quick-insert for database (vector of all generations to consider)
#' @param cohorts Quick-insert for database (vector of names of cohorts to consider)
#' @examples
#' data(ex_pop)
#' effects <- get.qtl.variance(ex_pop)
#' @return matrix with SNP / Chr / estimated effect variance
#' @export

get.qtl.variance <- function(population, gen=NULL, database=NULL, cohorts=NULL){
  if(length(gen)==0 && length(database)==0 && length(cohorts)==0){
    gen <- 1
  }
  geno <- get.geno(population, gen=gen, database=database, cohorts=cohorts)

  marker_variance <- list()
  for(index in 1:(length(population$info$real.bv.add)-1)){

    temp1 <- population$info$real.bv.add[[index]]
    if(length(temp1)>1){
      snp <- (temp1[,1] + c(0,population$info$cumsnp)[temp1[,2]])
      freq0 <- rowMeans(geno[snp,]==0)
      freq1 <- rowMeans(geno[snp,]==1)
      freq2 <- 1 - freq0 - freq1

      expect <- temp1[,3]*freq0 + temp1[,4] * freq1 + temp1[,5]* freq2
      expect2 <- temp1[,3]^2*freq0 + temp1[,4]^2 * freq1 + temp1[,5]^2* freq2

      marker_variance[[index]] <- cbind(temp1[,1:2], expect2 - expect^2)
    }


  }

  return(marker_variance)

}
