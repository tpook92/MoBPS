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

#' Devolopment of genetic/breeding value
#'
#' Function to plot genetic/breeding values for multiple generation/cohorts
#' @param population population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param json If TRUE extract which cohorts to plot according to the json-file used in json.simulation
#' @param ibd.obs Number of Individual pairs to sample for IBD estimation
#' @param hbd.obs Number of Individuals to sample for HBD estimation
#' @export

kinship.development <- function(population, database=NULL, gen=NULL, cohorts=NULL,
                           json=FALSE, ibd.obs=50, hbd.obs=10){

  if(json){
    ids <- to_plot <- numeric(length(population$info$json[[1]]))
    for(index in 1:length(ids)){
      ids[index] <- population$info$json[[1]][[index]]$label
      if(length(population$info$json[[1]][[index]]$'BV Plot')>0){
        to_plot[index] <- population$info$json[[1]][[index]]$'BV Plot'
      }
      cohorts <- ids[which(to_plot=="Yes")]
      if(length(cohorts)==0){
        cohorts <- population$info$cohorts[,1]
      }
    }
  }
  inbred <- NULL
  if(length(gen)>0){
    for(index in gen){
      inbred <- rbind(inbred, kinship.emp.fast(population=population, gen=index, ibd.obs=ibd.obs, hbd.obs=hbd.obs))
    }
  }
  if(length(database)>0){
    for(index in 1:nrow(database)){
      inbred <- rbind(inbred, kinship.emp.fast(population=population, database=database[index,1:2, drop=FALSE], ibd.obs=ibd.obs, hbd.obs=hbd.obs))
    }
  }
  if(length(cohorts)>0){
    for(index in cohorts){
      inbred <- rbind(inbred, kinship.emp.fast(population=population, cohorts=index, ibd.obs=ibd.obs, hbd.obs=hbd.obs))
    }
  }


  graphics::plot(inbred[,1], main=paste("Kinship development"), xaxt='n', xlab="", ylim=c(0,1), col="red", lwd=2, type="l", ylab="")
  graphics::lines(inbred[,2], main=paste("Kinship development"), xaxt='n', xlab="", ylim=c(0,1), lwd=2, col="blue")

  display.names <- NULL
  if(length(gen)>0){
    display.names <- c(display.names, paste("Generation", gen))
  }
  if(length(database)>0){
    display.names <- c(display.names, paste0(c("M", "F")[database[,2]], "_", database[,1]))
  }
  if(length(cohorts)>0){
    display.names <- c(display.names, cohorts)
  }
  graphics::axis(1, at=1:length(display.names), labels=display.names,las=2)
  graphics::legend("topleft",c("IBD", "HBD"), lty=c(1,1), col=c("red","blue"), lwd=c(2,2))

  return(inbred)
}





