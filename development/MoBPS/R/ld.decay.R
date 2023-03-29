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

#' Generate LD plot
#'
#' Generate LD pot
#' @param population Population list
#' @param genotype.dataset Genotype dataset (default: NULL - just to save computation time when get.geno was already run)
#' @param dist Manuel input of marker distances to analyse
#' @param step Stepsize to calculate LD
#' @param max Maximum distance between markers to consider for LD-plot
#' @param max.cases Maximum number of marker pairs to consider of each distance (default: 100; randomly sampled!)
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param chromosome Only consider a specific chromosome in calculations (default: 1)
#' @param type Compute LD decay according to following distance measure between markers (default: "snp", alt: "bp", "cM")
#' @param plot Set to FALSE to not generate an LD plot
#' @param xlim Axis limits for the x-axis in the LD plot
#' @param ylim Axis limits for the y-axis in the LD plot
#' @examples
#' data(ex_pop)
#' ld.decay(population=ex_pop, gen=5)
#' @return LD-decay plot for in gen/database/cohorts selected individuals
#' @export

ld.decay <- function(population, genotype.dataset=NULL, chromosome=1, dist =NULL,  step=5, max=500, max.cases = 100, database=NULL, gen=NULL, cohorts= NULL,
                     type="snp", plot=FALSE, xlim = NULL, ylim = NULL){
  max <- min(population$info$snp[chromosome]-1, max)
  if(length(genotype.dataset)==0){
    dataset <- t(get.geno(population, chromosome = chromosome, gen=gen, database=database,cohorts=cohorts))
  } else{
    dataset <- t(genotype.dataset)
  }

  if(length(dist)>0){
    calc <- dist
  } else{
    calc <- unique(c(1:(step-1),seq(from=step, to=max, by=step)))
  }
  if(type=="snp"){
    ld <- numeric(length(calc))
  } else{
    ld <- list()
    dlist <- list()
  }

  for(index in 1:length(calc)){


    cases <- 1:(population$info$snp[chromosome]-calc[index])
    if(length(cases)> max.cases){
      cases <- sample(cases, max.cases)
    }
    lds <- numeric(length(cases))
    temp1 <- 1
    for(index2 in cases){
      suppressWarnings(lds[temp1] <- stats::cor(dataset[,index2], dataset[,index2+calc[index]]))
      temp1 <- temp1 +1
    }
    if(type=="snp"){
      ld[index] <- mean(lds^2, na.rm=TRUE)
    } else{
      ld[[index]] <- lds^2
      if(type=="bp"){
        dlist[[index]] <- population$info$bp[cases + calc[index]] - population$info$bp[cases]
      } else if(type=="cM" || type=="cm"){
        dlist[[index]] <- population$info$snp.position[cases + calc[index]] - population$info$snp.position[cases]
      }

    }

  }

  if(type=="snp"){

    smooth1 <- stats::ksmooth(calc[1:(length(calc)/5)], ld[1:(length(calc)/5)], bandwidth = step*3, kernel="normal", x.points = calc[1:(length(calc)/5)])
    smooth2 <- stats::ksmooth(calc[-(1:(length(calc)/5))], ld[-(1:(length(calc)/5))], bandwidth = step*10 , kernel="normal", x.points = calc[-(1:(length(calc)/5))])
    smooth1$x[1] <- calc[1]
    smooth1$y[1] <- ld[1]
    smooth1$x <- c(smooth1$x, smooth2$x)
    smooth1$y <-c(smooth1$y, smooth2$y)

    if(plot){

      tryCatch(  {
        graphics::plot(calc, ld, xlab="distance in SNP", ylab=expression(r^2), main=paste0("LD structure on chromosome ", chromosome),
                       ylim = NULL, xlim = NULL)
        graphics::lines(smooth1, col="red", lwd=2)
      },
      error = function(e) {})


    }


    list(calc, ld, smooth1)
  } else{
    a <- unlist(ld)
    b <- unlist(dlist)

    b <- b[!is.na(a)]
    a <- a[!is.na(a)]
    evs <- numeric(length(dlist))
    for(index in 1:length(evs)){
      evs[index] <- mean(dlist[[index]])
    }
    smooth1 <- stats::ksmooth(b,a, x.points = evs, bandwidth = mean(diff(evs))*10, kernel="normal")
    smooth2 <- stats::ksmooth(b,a, x.points = evs, bandwidth = mean(diff(evs))*3, kernel="normal")
    smooth1$x[1:(length(calc)/5)] <- smooth2$x[1:(length(calc)/5)]
    if(type=="cm"||type=="cM"){
      type <- "Morgan"
    }
    if(plot){

      tryCatch(  {
        graphics::plot(smooth1 , xlab=paste0("distance in ", type), ylab=expression(r^2), main=paste0("LD structure on chromosome ", chromosome))
      },
      error = function(e) {})


    }


    list(a,b,smooth1)
  }



}
