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

#' Export underlying true breeding values
#'
#' Function to export underlying true breeding values
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param ids AA
#' @param import.position.calculation Function to calculate recombination point into adjacent/following SNP
#' @param decodeOriginsU Used function for the decoding of genetic origins [[5]]/[[6]]
#' @param plot Set TRUE to generate a visualization of genetic origins
#' @examples
#' data(ex_pop)
#' get.bv(ex_pop, gen=2)
#' @return Genomic value of in gen/database/cohorts selected individuals
#' @export


get.pool <- function(population, gen= NULL, database = NULL, cohorts = NULL, ids = NULL,
                      plot = FALSE, import.position.calculation=NULL, decodeOriginsU=decodeOriginsR){

  database <- get.database(population, gen = gen, database = database,
                           cohorts = cohorts, id = ids)

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
    gen = database[index5,1]
    sex = database[index5,2]
    for(nr in database[index5,3]:database[index5,4]){


      if(!is.na(population$breeding[[gen]][[sex+36]][[nr]])){
        segments[, cindex + 1:2] = population$breeding[[gen]][[sex+36]][[nr]]
      } else{
        for(hapnr in 1:2){
          current.recombi <- population$breeding[[gen]][[sex]][[nr]][[hapnr]]
          current.ursprung <- population$breeding[[gen]][[sex]][[nr]][[hapnr+4]]

          pool = numeric(length(current.ursprung))

          ursprung_type = unique(current.ursprung)
          for(index in 1:length(ursprung_type)){
            temp1 = decodeOriginsU(ursprung_type,index)
            pool[current.ursprung==ursprung_type[index]] = population$breeding[[temp1[1]]][[temp1[2]+36]][[temp1[3]]]
          }
          keep = c(TRUE, diff(pool)!=0)
          current.ursprung = current.ursprung[keep]
          current.recombi = current.recombi[c( keep, TRUE)]

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
            ursprung <- decodeOriginsU(ursprung_type,index2)
            ursprung[1] <- origin[ursprung[1]]

            relevant.snp.list = list()
            to_check = which(current.ursprung==ursprung_type[index2])
            for(index3 in 1:length(to_check)){
              tmp1 = to_check[index3]
              if((till[tmp1+1] >=(till[tmp1]+1))){
                relevant.snp.list[[index3]] <- (till[tmp1]+1):till[tmp1+1]
              }
            }
            relevant.snp = unlist(relevant.snp.list)

            segments[relevant.snp, hapnr + cindex] <-population$breeding[[ursprung[1]]][[ursprung[2]+36]][[ursprung[3]]]

          }

        }
      }

      cindex = cindex + 2

    }
  }


  if(plot){
    plot(0,-10, ylim=c(0,ncol(segments)), xlim=c(0, nrow(segments)), xlab="SNPs", ylab="haplotypes", yaxt="n")
    if(ncol(segments)==2){
      lab = c("Haplotype 1", "Haplotype 2")
    } else{
      lab = paste0("I", sort(rep(1:(ncol(segments)/2), 2)) , "H", 1:2)
    }
    graphics::axis(2, at=0.5+ seq(0, ncol(segments)-1, by=1), labels = lab)
    for(ind in 1:(ncol(segments)/2)){
      for(pair in 1:2){
        switch = c(1,which(c(0,diff(segments[,pair + ind*2-2]))!=0), nrow(segments)+1)
        for(index1 in 1:(length(switch)-1)){
          graphics::polygon(c(switch[index1], switch[index1+1], switch[index1+1], switch[index1]), c(pair-1+ ind*2-2, pair-1+ ind*2-2, pair+ ind*2-2, pair+ ind*2-2),
                            lty=0, col=segments[switch[index1],pair+ ind*2-2])
        }
      }
    }
    if(population$info$chromosome>1){
      graphics::abline(v=population$info$cumsnp[-length(population$info$cumsnp)], lty=3, lwd = 3, col = "grey")
    }


  }

  return(segments)
}
