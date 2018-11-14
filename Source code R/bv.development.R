'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2017 -- 2018  Torsten Pook

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

#' Entwicklung der verschiedenen Zuchtwerte
#'
#' Plot zur Entwicklung der realen/beobachteten/geschaetzten Zuchtwerte
#' @param population population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param confidence Draw confidence intervals for (1- bv, 2- bve, 3- pheno; default: c(1,2,3))
#' @param quantile Quantile of the confidence interval to draw (default: 0.05)
#' @param bvrow Which traits to display (for multiple traits separte plots (par(mfrow)))
#' @param ignore.zero Cohorts with only 0 individuals are not displayed (default: TRUE)
#' @export

bv.development <- function(population, database=NULL, gen=NULL, cohorts=NULL,
                           confidence= c(1,2,3), quantile=0.95, bvrow="all", ignore.zero=TRUE){

  if(length(bvrow)==1 && bvrow=="all"){
    bvrow <- 1:population$info$bv.nr
  }
  bv <- list()
  bve <- list()
  pheno <- list()

  graphics::par(mfrow=c(1,length(bvrow)))

  if(length(gen)>0){
    for(index in gen){
      bv[[length(bv)+1]] <- get.bv(population, gen=index)
      bve[[length(bve)+1]] <- get.bve(population, gen=index)
      pheno[[length(pheno)+1]] <- get.pheno(population, gen=index)
    }
  }
  if(length(database)>0){
    for(index in 1:nrow(database)){
      bv[[length(bv)+1]] <- get.bv(population, database=database[index,,drop=FALSE])
      bve[[length(bve)+1]] <- get.bve(population, database=database[index,,drop=FALSE])
      pheno[[length(pheno)+1]] <- get.pheno(population, database=database[index,,drop=FALSE])
    }
  }
  if(length(cohorts)>0){
    for(index in cohorts){
      bv[[length(bv)+1]] <- get.bv(population, cohorts=index)
      bve[[length(bve)+1]] <- get.bve(population, cohorts=index)
      pheno[[length(pheno)+1]] <- get.pheno(population, cohorts=index)
    }
  }


  for(nr in bvrow){
    color <- c("black", "blue", "red")
    q <- stats::qnorm(1 - (1-quantile)/2)
    means <- sds <- all0 <- matrix(NA, nrow=3, ncol=length(bv))

    for(index in 1:length(bv)){
      means[1,index] <- base::mean(bv[[index]][nr,])
      means[2,index] <- base::mean(bve[[index]][nr,])
      means[3,index] <- base::mean(pheno[[index]][nr,])
      sds[1,index] <- stats::sd(bv[[index]][nr,])
      sds[2,index] <- stats::sd(bve[[index]][nr,])
      sds[3,index] <- stats::sd(pheno[[index]][nr,])
      all0[1,index] <- base::prod(bv[[index]][nr,]==0)
      all0[2,index] <- base::prod(bve[[index]][nr,]==0)
      all0[3,index] <- base::prod(pheno[[index]][nr,]==0)
    }
    if(ignore.zero==TRUE){
      means[all0*(1:length(means))] <- NA
      sds[all0*(1:length(means))] <- NA
    }

    ub <- means - sds * q
    ob <- means + sds * q


    inc <- (rowSums(all0)!=ncol(all0))* (1:3)
    graphics::plot(means[1,],type="l",  main=paste("Development of breeding values - Trait", nr), xlab="groups", ylab="breeding value", ylim=c(min(ub, na.rm=TRUE),max(ob, na.rm=TRUE)), lwd=2)
    graphics::lines(means[2,], col="blue", lwd=2)
    graphics::lines(means[3,], col="red", lwd=2)
    for(art in confidence){
      graphics::lines(ob[art,], col=color[art], lty=3, lwd=2)
      graphics::lines(ub[art,], col=color[art], lty=3, lwd=2)
    }
    graphics::legend("bottomright",c("breeding value","phenotype","breeding value estimate")[inc], lty=c(1,1,1)[inc], col=color[inc], lwd=c(2,2,2)[inc])
  }


}
