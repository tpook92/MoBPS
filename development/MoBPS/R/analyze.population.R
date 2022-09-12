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

#' Analyze allele frequency of a single marker
#'
#' Analyze allele frequency of a single marker
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param chromosome Number of the chromosome of the relevant SNP
#' @param snp Number of the relevant SNP
#' @param snp.name Name of the SNP to analyze
#' @examples
#' data(ex_pop)
#' analyze.population(ex_pop, snp=1, chromosome=1, gen=1:5)
#' @return Frequency of AA/AB/BB in selected gen/database/cohorts
#' @export

analyze.population <- function(population, chromosome = NULL, snp=NULL, snp.name=NULL, database=NULL, gen=NULL, cohorts=NULL){

  if(length(snp.name)==1){
    n.snp <- which(population$info$snp.name == snp.name)
  } else{
    p.snp <- sum(population$info$length[0:(chromosome-1)]) + population$info$position[[chromosome]][snp]
    n.snp <- sum(population$info$snp[0:(chromosome-1)]) + snp
  }
  groups <- sum(nrow(database), length(gen), length(cohorts))
  state <- matrix(0, nrow=3, ncol = groups)


  col <- 1
  if(length(gen)>0){
    for(index in 1:length(gen)){
      genos <- get.geno(population, gen = gen[index])[n.snp,]
      state[,col] <- c(sum(genos==0), sum(genos==1), sum(genos==2))
      col <- col + 1
    }
  }
  if(length(database)>0){
    for(index in 1:nrow(database)){
      genos <- get.geno(population, database = database[index,,drop=FALSE])[n.snp,]
      state[,col] <- c(sum(genos==0), sum(genos==1), sum(genos==2))
      col <- col + 1
    }
  }
  if(length(cohorts)>0){
    for(index in 1:length(cohorts)){
      genos <- get.geno(population, cohorts = cohorts[index])[n.snp,]
      state[,col] <- c(sum(genos==0), sum(genos==1), sum(genos==2))
      col <- col + 1
    }
  }
  datatime <- c(gen, database[,1], as.numeric(population$info$cohorts[cohorts,2]))

  state.prob <- t(t(state)/colSums(state))
  maxp <- max(state.prob)
  graphics::plot(datatime ,state.prob[1,],xlim=c(min(datatime),max(datatime)),ylim=c(0,maxp),type="l",xlab="generation", ylab="frequency", lwd=3,
                 main="")
  graphics::lines(datatime ,state.prob[2,],lty=2, lwd=3)
  graphics::lines(datatime ,state.prob[3,],lty=3, lwd=3)
  graphics::points(datatime ,state.prob[1,], pch = 0)
  graphics::points(datatime ,state.prob[2,], pch = 1)
  graphics::points(datatime ,state.prob[3,], pch = 2)
  graphics::legend("topleft",legend = c("AA","AB","BB"),lty=c(1,2,3), lwd=c(3,3,3), pch=c(0,1,2))
  return(state)
}



