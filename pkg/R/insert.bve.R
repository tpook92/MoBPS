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

#' Export estimated breeding values
#'
#' Function to export estimated breeding values
#' @param population Population list
#' @param bves Matrix of breeding values to enter (one row per individual with 1 element coding individual name)
#' @param type which time of values to input (default: "bve", alt: "bv", "pheno")
#' @param na.override Set to TRUE to also enter NA values (Default: FALSE - those entries will be skipped)
#' @param count Counting for economic cost calculation (default: 1 - (one observation (for "pheno"), one genotyping (for "bve")))
#' @examples
#' data(ex_pop)
#' bv <- get.bv(ex_pop, gen=2)
#' new.bve <- cbind( colnames(bv), bv[,1]) ## Unrealistic but you do not get better than this!
#' ex_pop <- insert.bve(ex_pop, bves=new.bve)
#' @return Population-List with newly entered estimated breeding values
#' @export

insert.bve <- function(population, bves, type="bve", na.override = FALSE,  count=1){

  add <- 2
  if(type=="bv"){
    add <- 6
  } else if(type=="pheno"){
    add <- 8
  }

  if((ncol(bves)-1)!=population$info$bv.nr){
    stop("Number of traits entered does not match with population! \n Enter NA colums if you dont want to overwrite a trait")
  }

  for(index in 1:nrow(bves)){
    sex <- as.numeric(substr(bves[index,1], start=1, stop=1)=="F") + 1
    split <- strsplit(bves[index,1], split=c("_"))
    nr <- as.numeric(substr(split[[1]][1], start=2, stop=nchar(split[[1]][1])))
    gen <- as.numeric(split[[1]][2])
    if(na.override & add == 8){
      population$breeding[[gen]][[sex+add]][,nr] <- as.numeric(bves[index,-1])
    } else{
      population$breeding[[gen]][[sex+add]][,nr][!is.na(as.numeric(bves[index,-1]))] <- as.numeric(bves[index,-1])[!is.na(as.numeric(bves[index,-1]))]
    }

    if(add==2){
      population$breeding[[gen]][[sex]][[nr]][[16]] <- count
    } else if(add==8){
      if(count > 0 ){
        population$breeding[[gen]][[sex]][[nr]][[15]] <- (!is.na(population$breeding[[gen]][[sex+add]][,nr]))* count
      }
    }

  }

  return(population)
}
