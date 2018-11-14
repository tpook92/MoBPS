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


#' THIS NEEDs SOME REWORK BEFORE USAGE
#'
#' THIS NEEDs SOME REWORK BEFORE USAGE
#' @param population Bisherige Populationsliste
#' @param gen Generation in der Tiere hinzugefuegt werden sollen
#' @param sex Vektor mit deterministischen Angabe des Geschlechts (1-m, 2-w) jeden Tieres (default: NULL)
#' @param sex.odds Wahrscheinlichkeit eines weiblichen Tieres
#' @param dataset Datensatz neuer Tiere (Plinkfile - Haplotypen)
#' @param new.class Migrationslevel der neuen Tiere
#' @param first.col Spalte im Datensatz mit dem ersten Tier
#' @param sex.s Gender of newly added individuals
#' @export

creation.add <- function(population, gen=1, sex.s=NULL, sex.odds=0.5, dataset, new.class = 0 , first.col=1 ){
  if(length(population$breeding[[1]][[1]])>0){
    dummy.animal <- population$breeding[[1]][[1]][[1]]
  } else{
    dummy.animal <- population$breeding[[1]][[2]][[1]]
  }
  if (requireNamespace("miraculix", quietly = TRUE)) {
    codeOriginsU <- miraculix::codeOrigins
    decodeOriginsU <- miraculix::decodeOrigins
  } else{
    codeOriginsU <- codeOriginsR
    decodeOriginsU <- decodeOriginsR
  }

  # compute origin.gen
  if(length(population$info$origin.gen)>0){
    take <- which(population$info$origin.gen==gen)
    if(length(take)==1){
      origin_code <- population$info$origin.gen[take]
    } else{
      if(length(population$info$origin.gen)<64){
        population$info$origin.gen <- c(population$info$origin.gen,as.integer(gen))
        origin_code <- length(population$info$origin.gen)
      } else{
        print("To many origin generation!")
        print("Delete second lowest origin.gen")
        switch <- sort(population$info$origin.gen, index.return=TRUE)[[2]]
        population$info$origin.gen[switch] <- as.integer(gen)
        origin_code <- switch
      }
    }
  } else{
    origin_code <- gen
  }
  new.animals <- (ncol(dataset)-first.col+1)/2
  if(length(sex.s)==0){
    sex.s <- stats::rbinom(new.animals,1, sex.odds)+1
  } else{
    sex.s <- sex.s
  }
  breeding.size <- c(sum(sex.s==1), sum(sex.s==2))
  if(length(population$breeding)<gen){
    population$breeding[[gen]] <- list()
    population$info$size <- rbind(population$info$size,0)
  }
  current.size <- c(1,1)
  for(sex in 1:2){
    if(length(population$breeding[[gen]])<=2 || length(population$breeding[[gen]][[sex]])==0){
      population$breeding[[gen]][[2+sex]] <- matrix(0, nrow=population$info$bv.nr, ncol=breeding.size[sex])
      population$breeding[[gen]][[4+sex]] <- rep(new.class, breeding.size[sex])
      population$breeding[[gen]][[6+sex]] <- matrix(0, nrow=population$info$bv.nr, ncol=breeding.size[sex])
      population$breeding[[gen]][[8+sex]] <- matrix(0, nrow=population$info$bv.nr, ncol=breeding.size[sex])
      #    } else if(length(population$breeding[[gen]][[sex+2]])==0){
      #      population$breeding[[gen]][[2+sex]] <- rep(0, breeding.size[sex])
      #      population$breeding[[gen]][[4+sex]] <- rep(new.class, breeding.size[sex])
      #     population$breeding[[gen]][[6+sex]] <- new.bv[sex,]
      #      population$breeding[[gen]][[8+sex]] <- new.bv.approx[sex,]
    } else{
      current.size[sex] <- length(population$breeding[[gen]][[4+sex]]) + 1
      population$breeding[[gen]][[2+sex]] <- cbind(population$breeding[[gen]][[sex+2]], matrix(0, nrow= population$info$bv.nr, ncol=breeding.size[sex]))
      population$breeding[[gen]][[4+sex]] <- c(population$breeding[[gen]][[sex+4]], rep(new.class, breeding.size[sex]))
      population$breeding[[gen]][[6+sex]] <- cbind(population$breeding[[gen]][[6+sex]], matrix(0, nrow= population$info$bv.nr, ncol=breeding.size[sex]))
      population$breeding[[gen]][[8+sex]] <- cbind(population$breeding[[gen]][[8+sex]], matrix(0, nrow= population$info$bv.nr, ncol=breeding.size[sex]))
    }
    if(length(population$breeding[[gen]][[sex]])==0){
      population$breeding[[gen]][[sex]] <- list()
    }
  }

  current.size <- c(length(population$breeding[[gen]][[1]]),length(population$breeding[[gen]][[2]]))+1

  current.row <- first.col
  for(index in 1:new.animals){
    sex <- sex.s[index]
    population$breeding[[gen]][[sex]][[current.size[sex]]] <- dummy.animal
    population$breeding[[gen]][[sex]][[current.size[sex]]][[5]] <- codeOriginsU(matrix(c(origin_code, sex,  current.size[sex], 1),nrow=(length(population$breeding[[gen]][[sex]][[current.size[sex]]][[1]])-1), ncol=4, byrow=TRUE))
    population$breeding[[gen]][[sex]][[current.size[sex]]][[6]] <- codeOriginsU(matrix(c(origin_code, sex,  current.size[sex], 2),nrow=(length(population$breeding[[gen]][[sex]][[current.size[sex]]][[2]])-1), ncol=4, byrow=TRUE))
    population$breeding[[gen]][[sex]][[current.size[sex]]][[7]] <- c(gen, sex, current.size[sex])
    population$breeding[[gen]][[sex]][[current.size[sex]]][[8]] <- c(gen, sex, current.size[sex])
    population$breeding[[gen]][[sex]][[current.size[sex]]][[9]] <- 1 - (dataset[,current.row]==population$info$snp.base[1,])
    population$breeding[[gen]][[sex]][[current.size[sex]]][[10]] <- 1 - (dataset[,(current.row+1)]==population$info$snp.base[1,])
    current.row <- current.row +2
    current.size[sex] <- current.size[sex] + 1
    population$info$size[gen,sex] <- population$info$size[gen,sex] +1
  }

  return(population)
}
