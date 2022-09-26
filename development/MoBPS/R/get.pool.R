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
#' @param plot Set TRUE to generate a visualization of genetic origins
#' @examples
#' data(ex_pop)
#' get.bv(ex_pop, gen=2)
#' @return Genomic value of in gen/database/cohorts selected individuals
#' @export


get.pool <- function(population, gen= NULL, database = NULL, cohorts = NULL, ids = NULL,
                     plot = FALSE){
  database <- get.database(population, gen = gen, database = database,
                           cohorts = cohorts, ids = ids)

  recombi = get.recombi(population, database = database)

  segments = matrix(nrow=sum(population$info$snp), ncol=length(recombi)*2)

  for(index in 1:length(recombi)){

    for(pair in 1:2){
      temp1 <- recombi[[index]][[pair+2]]
      for(rec in 1:(length(recombi[[index]][[pair]])-1)){
        activ_rec = (population$info$snp.position < recombi[[index]][[pair]][rec+1]) & (population$info$snp.position > recombi[[index]][[pair]][rec])
        pool = population$breeding[[temp1[rec,1]]][[temp1[rec,2]+36]][[temp1[rec,3]]]
        segments[activ_rec,index*2+pair-2] <- pool
      }

    }


  }

  if(plot){
    plot(0,-10, ylim=c(0,ncol(segments)), xlim=c(0, nrow(segments)), xlab="SNPs", ylab="haplotypes", yaxt="n")
    if(ncol(segments)==2){
      lab = c("Haplotype 1", "Haplotype 2")
    } else{
      lab = paste0("I", sort(rep(1:(ncol(segments)/2), 2)) , "H", 1:2)
    }
    axis(2, at=0.5+ seq(0, ncol(segments)-1, by=1), labels = lab)
    for(ind in 1:(ncol(segments)/2)){
      for(pair in 1:2){
        switch = c(1,which(c(0,diff(segments[,pair + ind*2-2]))!=0), nrow(segments)+1)
        for(index1 in 1:(length(switch)-1)){
          polygon(c(switch[index1], switch[index1+1], switch[index1+1], switch[index1]), c(pair-1+ ind*2-2, pair-1+ ind*2-2, pair+ ind*2-2, pair+ ind*2-2),
                  lty=0, col=segments[switch[index1],pair+ ind*2-2])
        }
      }
    }

  }


}

polygon()
