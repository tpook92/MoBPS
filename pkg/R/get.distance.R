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

#' Calculate Nei distance between two or more population
#'
#' Function to calculate Nei's distance between two or more population
#' @param population population list
#' @param type Chose type of distance to compute (default: Neis standard genetic distance "nei"). Alt: Reynolds distance ("reynold"), Cavalli-Sforza ("cavalli"), Neis distance ("nei_distance"), Neis minimum distance ("nei_minimum")
#' @param marker Vector with SNPs to consider (Default: "all" - use of all markers)
#' @param per.marker Set to TRUE to return per marker statistics on genetic distances
#' @param database1 First Groups of individuals to consider
#' @param gen1 Quick-insert for database (vector of all generations to consider)
#' @param cohorts1 Quick-insert for database (vector of names of cohorts to consider)
#' @param database2 Second Groups of individuals to consider
#' @param gen2 Quick-insert for database (vector of all generations to consider)
#' @param cohorts2 Quick-insert for database (vector of names of cohorts to consider)
#' @param database.list List of databases to consider (use when working with more than 2 populations)
#' @param gen.list Quick-insert for database (vector of all generations to consider)
#' @param cohorts.list Quick-insert for database (vector of names of cohorts to consider)
#'
#' @examples
#' data(ex_pop)
#' get.distance(ex_pop, database1 = cbind(1,1), database2 = cbind(1,2))
#' @return Population list
#' @export

get.distance <- function(population, type="nei", marker = "all", per.marker=FALSE,
                              gen1 = NULL, database1 = NULL, cohorts1 = NULL,
                             gen2 = NULL, database2 = NULL, cohorts2= NULL,
                             database.list = NULL, gen.list = NULL, cohorts.list = NULL){

  if(length(marker)==1 && marker=="all"){
    marker <- 1:sum(population$info$snp)
  }
  database1 <- get.database(population, gen1, database1, cohorts1)
  database2 <- get.database(population, gen2, database2, cohorts2)

  if(!is.list(database.list)){
    database.list <- list()
  }
  if(!is.list(gen.list)){
    gen.list <- list()
  }
  if(!is.list(cohorts.list)){
    cohorts.list <- list()
  }

  if(length(database1)>0){
    database.list[[1]] <- database1
  }
  if(length(database2)>0){
    database.list[[2]] <- database2
  }

  for(index in 1:max(length(database.list), length(gen.list), length(cohorts.list))){
    database.list[[index]] <- get.database(population,
                                           if(length(gen.list)>=index){gen.list[[index]]},
                                           if(length(database.list)>=index){database.list[[index]]},
                                           if(length(cohorts.list)>=index){cohorts.list[[index]]})
  }

  p_i <- list()

  n_p <- length(database.list)

  for(index in 1:n_p){
    p_i[[index]] <- rowMeans(get.geno(population, database = database.list[[index]])) / 2
  }



  D <- matrix(0, ncol=n_p, nrow=n_p)
  if(n_p>2 && per.marker){
    per.marker <- FALSE
    warning("Per marker statistic only available for exactly two populations")
  }

  for(index in 1:(n_p-1)){
    for(index2 in (index+1):n_p){
      p1 <- p_i[[index]][marker]
      p2 <- p_i[[index2]][marker]
      temp1 <- NULL
      if(type=="nei"){
        temp1 <-  -log( ( (p1*p2 + (1-p1)*(1-p2)) ) / sqrt((p1^2 + (1-p1)^2) * (p2^2+(1-p2)^2))     )
      } else if(type=="cavalli"){
        temp1 <-  2 / pi * sqrt(2 * (1- sqrt(p1*p2) - sqrt((1-p1)*(1-p2)) ))
      } else if(type=="reynold"){
        temp1 <-  (p1-p2)^2 / (1 - p1*p2 - (1-p1)*(1-p2))
      } else if(type=="nei_minimum"){
        temp1 <- (p1-p2)^2
      } else if(type=="nei_distance"){
        temp1 <- 1 - sqrt(p1*p2) - sqrt((1-p1) * (1-p2))
      }
        D[index,index2] <- D[index2,index] <-sum(temp1, na.rm=TRUE) / length(p1)
    }
  }


  if(per.marker){
    return(temp1)
  } else{
    return(D)
  }

}




