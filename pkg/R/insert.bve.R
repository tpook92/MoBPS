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

#' Manually enter estimated breeding values
#'
#' Function to manually enter estimated breeding values
#' @param population Population list
#' @param bves Matrix of breeding values to enter (one row per individual with 1 element coding individual name)
#' @param type which time of values to input (default: "bve", alt: "bv", "pheno")
#' @param na.override Set to TRUE to also enter NA values (Default: FALSE - those entries will be skipped)
#' @param count Counting for economic cost calculation (default: 1 - (one observation (for "pheno"), one genotyping (for "bve")))
#' @param count.only.increase Set to FALSE to reduce the number of observation for a phenotype to "count" (default: TRUE)
#' @param all.copys Set FALSE to only enter target values for a specific copy of the individual (default: TRUE)
#' @examples
#' data(ex_pop)
#' bv <- get.bv(ex_pop, gen=2)
#' new.bve <- cbind( colnames(bv), bv[,1]) ## Unrealistic but you do not get better than this!
#' ex_pop <- insert.bve(ex_pop, bves=new.bve)
#' @return Population-List with newly entered estimated breeding values
#' @export

insert.bve <- function(population, bves, type="bve", na.override = FALSE,  count=1, count.only.increase=TRUE, all.copys = TRUE){

  add <- 2
  if(type=="bv"){
    add <- 6
  } else if(type=="pheno"){
    add <- 8
  } else if(type=="reli"){
    add <- 16
  }

  # check if the old formatting is used:

  suppressWarnings({
    old_format = (ncol(bves)-1)==population$info$bv.nr && ((substr(bves[1,1], start=1, stop=1) %in% c("F", "M")) || (!is.na(as.numeric(bves[1,1]) && (as.integer(bves[1,1]) == as.numeric(bves[1,1])))))
  })


  if(!old_format){
    if(nrow(bves) == population$info$bv.nr && (ncol(bves)-1)!=population$info$bv.nr){
      bves = cbind(colnames(bves), t(bves))
    }

    if(ncol(bves) == population$info$bv.nr){
      bves = cbind(rownames(bves), bves)
    }
  }



  if((ncol(bves)-1)!=population$info$bv.nr){
    stop("Number of traits entered does not match with population! \n Enter NA colums if you dont want to overwrite a trait")
  }

  is_id <- !(substr(bves[1,1], start=1, stop=1) %in% c("F", "M"))

  if(is_id){

    position = get.database(population, id = bves[,1], keep.order = TRUE)[,1:3,drop=FALSE]
    sex = position[,2]
    gen = position[,1]
    nr = position[,3]

  } else{
    sex <- as.numeric(substr(bves[,1], start=1, stop=1)=="F")+1
    split <- unlist(strsplit(bves[,1], split=c("_")))
    gen <- as.numeric(split[1:length(sex)*2])
    nr <- as.numeric(substr(split[1:length(sex)*2-1], start=2, stop=100))
  }


  # check if the dataset in general contains duplicates
  if(all.copys && prod(1==as.numeric(population$info$cohorts[-1,10]) - as.numeric(population$info$cohorts[-nrow(population$info$cohorts),11]))!=1){

    sex_list = list()
    gen_list = list()
    nr_list = list()
    index_list = list()
    for(index in 1:length(nr)){
      list_copies = population$breeding[[gen[index]]][[sex[index]]][[nr[index]]][[21]]
      if(nrow(list_copies) >= 1){
        gen_list[[index]] = list_copies[,1]
        sex_list[[index]] = list_copies[,2]
        nr_list[[index]] = list_copies[,3]
        index_list[[index]] = rep(index, nrow(list_copies))
      }
    }
    sex = unlist(sex_list)
    gen = unlist(gen_list)
    nr = unlist(nr_list)
    bves = bves[unlist(index_list),]

  }


  databases <- cbind(gen, sex, nr)

  groups <- unique(databases[,1:2])

  for(index in 1:nrow(groups)){
    activ <- databases[,1]==groups[index,1] & databases[,2]==groups[index,2]
    database_activ <- databases[activ,,drop = FALSE]
    bves_activ <- t(bves[activ,-1, drop = FALSE])
    storage.mode(bves_activ) <- "numeric"

    if(sum(diff(database_activ[,3])==1)==(nrow(database_activ)-1)){
      from <- min(database_activ[,3])
      to <- max(database_activ[,3])

      if(na.override){
        population$breeding[[groups[index,1]]][[groups[index,2]+add]][,from:to] <- bves_activ
      } else{
        population$breeding[[groups[index,1]]][[groups[index,2]+add]][,from:to][!is.na(bves_activ)] <- bves_activ[!is.na(bves_activ)]
      }

    } else{

      if(na.override){
        for(index2 in 1:nrow(database_activ)){
          nr <- database_activ[index2,3]
          population$breeding[[groups[index,1]]][[groups[index,2]+add]][,nr] <- bves_activ[,index2]
        }
      } else{
        for(index2 in 1:nrow(database_activ)){
          nr <- database_activ[index2,3]
          population$breeding[[groups[index,1]]][[groups[index,2]+add]][,nr][!is.na(bves_activ[,index2]) ]  <- as.numeric(bves_activ[,index2])[!is.na(bves_activ[,index2]) ]
        }
      }


    }

    if(add==8){

      for(index2 in 1:nrow(database_activ)){
        sex <- groups[index,2]
        gen <- groups[index,1]
        nr <- database_activ[index2,3]

        temp1 <- (!is.na(population$breeding[[gen]][[sex+add]][,nr]))* count
        if(count.only.increase){
          population$breeding[[gen]][[sex]][[nr]][[15]][population$breeding[[gen]][[sex]][[nr]][[15]]<temp1] <- temp1[population$breeding[[gen]][[sex]][[nr]][[15]]<temp1]
        } else{
          population$breeding[[gen]][[sex]][[nr]][[15]] <- temp1
          if(length(population$breeding[[gen]][[sex]][[nr]][[24]])>0 && ncol(population$breeding[[gen]][[sex]][[nr]][[24]])>max(temp1)){
            if(max(temp1)>0){
              population$breeding[[gen]][[sex]][[nr]][[24]] <- population$breeding[[gen]][[sex]][[nr]][[24]][,1:max(temp1)]
            } else{
              population$breeding[[gen]][[sex]][[nr]][24] <- list(NULL) ## Only single bracket to not reduce length of list
            }

          }

        }

        for(bven in which(temp1==1)){
          if(length(population$breeding[[gen]][[sex]][[nr]][[27]]) < bven || length(population$breeding[[gen]][[sex]][[nr]][[27]][[bven]])<=1){
            population$breeding[[gen]][[sex]][[nr]][[27]][[bven]] = population$breeding[[gen]][[sex+add]][bven,nr]
          }
        }
      }

    }

  }

  return(population)
}

#' Manually enter phenotypes
#'
#' Function to manually enter phenotypes
#' @param population Population list
#' @param phenos Matrix of phenotypes to enter (one row per individual with 1 element coding individual name)
#' @param na.override Set to TRUE to also enter NA values (Default: FALSE - those entries will be skipped)
#' @param count Counting for economic cost calculation (default: 1 - (one observation (for "pheno"), one genotyping (for "bve")))
#' @param count.only.increase Set to FALSE to reduce the number of observation for a phenotype to "count" (default: TRUE)
#' @param all.copys Set FALSE to only enter target values for a specific copy of the individual (default: TRUE)
#' @examples
#' data(ex_pop)
#' bv <- get.bv(ex_pop, gen=2, use.id = FALSE)
#' new.bve <- cbind( colnames(bv), bv[,1]) ## Unrealistic but you do not get better than this!
#' ex_pop <- insert.pheno(ex_pop, phenos=new.bve)
#' @return Population-List with newly entered phenotypes
#' @export
#'
insert.pheno <- function(population, phenos, na.override = FALSE,  count=1, count.only.increase=TRUE, all.copys = TRUE){
  population <- insert.bve(population, bves=phenos, type="pheno", na.override=na.override, count=count,
                           count.only.increase = count.only.increase, all.copys = all.copys)
}

#' Manually enter breeding values
#'
#' Function to manually enter breeding values
#' @param population Population list
#' @param bvs Matrix of phenotypes to enter (one row per individual with 1 element coding individual name)
#' @param na.override Set to TRUE to also enter NA values (Default: FALSE - those entries will be skipped)
#' @param count Counting for economic cost calculation (default: 1 - (one observation (for "pheno"), one genotyping (for "bve")))
#' @param count.only.increase Set to FALSE to reduce the number of observation for a phenotype to "count" (default: TRUE)
#' @param all.copys Set FALSE to only enter target values for a specific copy of the individual (default: TRUE)
#' @examples
#' data(ex_pop)
#' bv <- get.bv(ex_pop, gen=2, use.id = FALSE)
#' new.bve <- cbind( colnames(bv), bv[,1]) ## Unrealistic but you do not get better than this!
#' ex_pop <- insert.bv(ex_pop, bvs=new.bve)
#' @return Population-List with newly entered breeding values
#' @export
#'
insert.bv <- function(population, bvs, na.override = FALSE,  count=1, count.only.increase=TRUE, all.copys = TRUE){
  population <- insert.bve(population, bves=bvs, type="bv", na.override=na.override, count=count,
                           count.only.increase = count.only.increase, all.copys = all.copys)
}
