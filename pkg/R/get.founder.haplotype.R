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

#' Derive founder haplotypes of selected individuals
#'
#' Function to derive founder haplotypes of selected individuals
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param id Individual IDs to search/collect in the database
#' @param plot Set to TRUE to visualize founder haplotypes (default: FALSE)
#' @param import.position.calculation XXX
#' @param decodeOriginsU XXX
#' @param use.id XXX
#' @param array XXX
#' @examples
#' data(ex_pop)
#' haplo <- get.haplo(ex_pop, gen=2)
#' @return Founder Haplotype data for in gen/database/cohorts selected individuals
#' @export

get.founder.haplotype <- function(population, gen= NULL, database = NULL, cohorts = NULL, id = NULL,
                     plot = FALSE, import.position.calculation=NULL, decodeOriginsU=decodeOriginsR,
                     use.id = TRUE, array = NULL){

  if(length(array)>0){
    temp1 <- which(population$info$array.name == array)
    if(length(temp1)==1){
      array <- temp1
    }
    if(!is.numeric(array)){
      stop("Input for array can not be assigned! Check your input!")
    }
  }

  database <- get.database(population, gen = gen, database = database,
                           cohorts = cohorts, id = id)

  nindi = sum(database[,4]-database[,3]+1)

  snp_temp = population$info$snp
  n.snps <- sum(snp_temp)
  equi = sum(population$info$snps.equidistant==TRUE)>0
  origin = population$info$origin.gen
  if(equi){
    lec = population$info$length.total
    n_chr = length(lec) -1
    between2_vec = length_chromo_vec = position_chromo_vec = numeric(n_chr)

    for(chr in 1:n_chr){
      between2_vec[chr] <- population$info$position[[chr]][1]
      length_chromo_vec[chr] <- sum(population$info$length[0:(chr-1)])
      position_chromo_vec[chr] <- sum(snp_temp[0:(chr-1)])
    }

  }
  last = population$info$snp.position[n.snps]
  segments = matrix(nrow=n.snps, ncol=nindi*2)



  cindex = 0
  for(index5 in 1:nrow(database)){
    for(nr in database[index5,3]:database[index5,4]){
      gen = database[index5,1]
      sex = database[index5,2]

      tmp = population$breeding[[gen]][[sex]][[nr]][[21]]
      if(length(tmp) > 3){
        gen = tmp[1]
        sex = tmp[2]
        nr = tmp[3]
      }


        for(hapnr in 1:2){
          current.recombi <- population$breeding[[gen]][[sex]][[nr]][[hapnr]] / population$info$recombination.precision
          current.ursprung <- population$breeding[[gen]][[sex]][[nr]][[hapnr+4]]

          keep = c(TRUE, diff(current.ursprung)!=0)
          current.ursprung = current.ursprung[keep]
          current.recombi = current.recombi[c( keep, TRUE)]

          pool = current.ursprung


          if(is.function(import.position.calculation)){
            till <- import.position.calculation(current.recombi)
          } else if(equi){
            till <- rep(0, length(current.recombi))
            for(index in 1:length(current.recombi)){
              current.chromo <- sum(current.recombi[index] >= lec)
              if(current.recombi[index]>last){
                till[index] <- n.snps
              } else{
                between2 <- between2_vec[current.chromo]
                length_chromo <- length_chromo_vec[current.chromo]
                position_chromo <- position_chromo_vec[current.chromo]
                till[index] <- ceiling((current.recombi[index] - length_chromo - between2) / (between2*2)) + position_chromo

              }
            }
          } else{
            till <- rep(0, length(current.recombi))
            for(index in 1:length(current.recombi)){
              till[index] <- find.snpbefore(current.recombi[index], population$info$snp.position)
              #        till[index] <- sum(population$info$snp.position <= current.recombi[index])
            }

          }

          ursprung_type = unique(current.ursprung) # due to clean up this needs to be recalculated

          for(index2 in 1:(length(ursprung_type))){
            ursprung <- ursprung_type[index2]
            relevant.snp.list = list()
            to_check = which(current.ursprung==ursprung_type[index2])
            for(index3 in 1:length(to_check)){
              tmp1 = to_check[index3]
              if((till[tmp1+1] >=(till[tmp1]+1))){
                relevant.snp.list[[index3]] <- (till[tmp1]+1):till[tmp1+1]
              }
            }
            relevant.snp = unlist(relevant.snp.list)

            segments[relevant.snp, hapnr + cindex] <- ursprung

          }

        }
      cindex = cindex + 2

      }

  }

  before = 0
  names <- numeric(nindi)
  for(row in 1:nrow(database)){
    animals <- database[row,]
    nanimals <- database[row,4] - database[row,3] +1
    if(nanimals>0){
      names[(before+1):(before+nanimals)] <- paste(if(animals[2]==1) "M" else "F", animals[3]:animals[4],"_", animals[1], sep="")
      before <- before + nanimals
    }
  }

  if(use.id){
    tmp1 <- get.id(population, database = database)
  } else{
    tmp1 <- names
  }

  tmp1 = paste0(rep(tmp1, each = 2), c("_1", "_2"))

  colnames(segments) = tmp1

  if(length(array)>0){
    segments = segments[population$info$array.markers[[array]],]

  }

  return(segments)
}









