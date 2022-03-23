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

#' Empirical kinship
#'
#' Function to compute empirical kinship for a set of individuals)
#' @param population Population list
#' @param animals List of animals to compute kinship for
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param sym If True derive matrix entries below principle-diagonal
#' @param use.id Set TRUE to use animal IDs for column/row-names in the output matrix
#' @examples
#' data(ex_pop)
#' kinship <- kinship.emp(population=ex_pop, database=cbind(2,1,1,25))
#' @return Empirical kinship matrix (IBD-based since Founders)
#' @export

kinship.emp <- function(population=NULL, animals=NULL, gen=NULL, database=NULL, cohorts=NULL, sym=FALSE, use.id=FALSE){


  if(length(animals)==0 || length(gen)>0 || length(database)>0 || length(cohorts)>0){
    database <- get.database(population, gen, database, cohorts)
    animals <- list()
    for(index in 1:nrow(database)){
      if(diff(database[index,3:4])>=0){
        for(index2 in database[index,3]:database[index,4]){
          animals[[length(animals)+1]] <- population$breeding[[database[index,1]]][[database[[index,2]]]][[index2]]
        }
      }

    }
  }
  n <- length(animals)
  kinship <- matrix(0, nrow=n, ncol=n)

  col_names <- get.id(population, database=database)
  if(use.id==FALSE){
    col_names <- names(col_names)
  }


  chrom.length <- max(animals[[1]][[1]])
  for( i in 1:n){
    for( j in i:n){
        chr <- list()
        chr[[1]] <- animals[[i]][[1]][-1]
        chr[[2]] <- animals[[i]][[2]][-1]
        chr[[3]] <- animals[[j]][[1]][-1]
        chr[[4]] <- animals[[j]][[2]][-1]
        origin <- list()
        origin[[1]] <- animals[[i]][[5]]
        origin[[2]] <- animals[[i]][[6]]
        origin[[3]] <- animals[[j]][[5]]
        origin[[4]] <- animals[[j]][[6]]

        activ <- c(1,1,1,1)
        prev <- 0
        activ.recom <- c(chr[[1]][activ[1]], chr[[2]][activ[2]], chr[[3]][activ[3]], chr[[4]][activ[4]])
        activ.ursprung <- c(origin[[1]][activ[1]],origin[[2]][activ[2]],origin[[3]][activ[3]],origin[[4]][activ[4]])



        for(steps in 1:(length(c(chr[[1]], chr[[2]], chr[[3]], chr[[4]]))-3)){
          activ.min <- which.min(activ.recom)[1]
          activ.posi <- chr[[activ.min]][activ[activ.min]]

          ibd.factor <- (sum(activ.ursprung[1] == activ.ursprung[3:4]) + sum(activ.ursprung[2] == activ.ursprung[3:4]))/4

          #ibd <- length(unique(activ.ursprung)) # Nur vergleich des Neuen mit bisherigen Rechenzeiteffizienter!
          #ibd.factor <- 1-ibd*0.25 + 0.25*(ibd==1)
          kinship[i,j] <- kinship[i,j] + ibd.factor * (activ.posi - prev) / chrom.length
          prev <- activ.posi

          activ[activ.min] <- min(activ[activ.min] +1, length(chr[[activ.min]]))
          activ.recom[activ.min] <- chr[[activ.min]][activ[activ.min]]
          activ.ursprung[activ.min] <-  origin[[activ.min]][activ[activ.min]]

        }
      }

  }
  if(sym==TRUE){
    for(i in 1:n){
      for(j in 1:i){
        kinship[i,j] <- kinship[j,i]
      }
    }
  }
  colnames(kinship) <- rownames(kinship) <- col_names

  return(kinship)
}

#' Approximate empirical kinship
#'
#' Function to compute empirical kinship for a set of individuals (not all pairs of individuals are evaluated)
#' @param population Population list
#' @param animals List of animals to compute kinship for
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param ibd.obs Number of Individual pairs to sample for IBD estimation
#' @param hbd.obs Number of Individuals to sample for HBD estimation
#' @examples
#' data(ex_pop)
#' kinship.emp.fast(population=ex_pop,gen=2)
#' @return Empirical kinship matrix (IBD-based since Founders) per gen/database/cohort
#' @export
#'
kinship.emp.fast <- function(population=NULL, animals=NULL, gen=NULL, database=NULL, cohorts=NULL, ibd.obs=50, hbd.obs=10){
  if(length(animals)==0 || length(gen)>0 || length(database)>0 || length(cohorts)>0){
    database <- get.database(population, gen=gen, database=database, cohorts=cohorts)
    animals <- list()
    for(index in 1:nrow(database)){
      if(diff(database[index,3:4])>=0){
        for(index2 in database[index,3]:database[index,4]){
          animals[[length(animals)+1]] <- population$breeding[[database[index,1]]][[database[[index,2]]]][[index2]]
        }
      }
    }
  }

  n <- length(animals)
  if(n==0){
    return(c(0,0))
  }

  chrom.length <- max(animals[[1]][[1]])

  if(chrom.length=="removed"){
    return(c(0,0))
  }

  if(n^2 <= (ibd.obs + hbd.obs)){
    i1 <- c(rep(1:n, n)[-(1:n + (1:n) * n - n )], 1:n)
    j1 <- c(sort(rep(1:n, n))[-(1:n + (1:n) * n - n )], 1:n)
    ibd.obs <- n * (n-1)
    hbd.obs <- n
  } else{
    i1 <- sample(1:n, ibd.obs+ hbd.obs, replace=if((ibd.obs+hbd.obs)>n){TRUE} else{FALSE})
    j1 <- sample(1:n, ibd.obs+ hbd.obs, replace=if((ibd.obs+hbd.obs)>n){TRUE} else{FALSE})
  }

  if(ibd.obs>0){
    while(sum(i1==j1)>0){
      j1[which(i1==j1)] <- sample(1:n, length(which(i1==j1)), replace=if((ibd.obs)>n){TRUE} else{FALSE})
    }
  }
  if(hbd.obs>0){
    j1[(ibd.obs+1):(ibd.obs+hbd.obs)] <- i1[(ibd.obs+1):(ibd.obs+hbd.obs)]
  }

  score <- numeric(length(i1))
  for(index in 1:length(i1)){
    i <- i1[index]
    j <- j1[index]
    chr <- list()
    chr[[1]] <- animals[[i]][[1]][-1]
    chr[[2]] <- animals[[i]][[2]][-1]
    chr[[3]] <- animals[[j]][[1]][-1]
    chr[[4]] <- animals[[j]][[2]][-1]
    origin <- list()
    origin[[1]] <- animals[[i]][[5]]
    origin[[2]] <- animals[[i]][[6]]
    origin[[3]] <- animals[[j]][[5]]
    origin[[4]] <- animals[[j]][[6]]

    activ <- c(1,1,1,1)
    prev <- 0
    activ.recom <- c(chr[[1]][activ[1]], chr[[2]][activ[2]], chr[[3]][activ[3]], chr[[4]][activ[4]])
    activ.ursprung <- c(origin[[1]][activ[1]],origin[[2]][activ[2]],origin[[3]][activ[3]],origin[[4]][activ[4]])



    for(steps in 1:(length(c(chr[[1]], chr[[2]], chr[[3]], chr[[4]]))-3)){
      activ.min <- which.min(activ.recom)[1]
      activ.posi <- chr[[activ.min]][activ[activ.min]]

      if((activ.posi - prev)>0){
        #ibd <- length(unique(activ.ursprung)) # Nur vergleich des Neuen mit bisherigen Rechenzeiteffizienter!
        #ibd <- sum(!duplicated(activ.ursprung)) # Nur vergleich des Neuen mit bisherigen Rechenzeiteffizienter!
        #ibd.factor <- 1-ibd*0.25 + 0.25*(ibd==1)

        ibd.factor <- (sum(activ.ursprung[1] == activ.ursprung[3:4]) + sum(activ.ursprung[2] == activ.ursprung[3:4]))/4
        score[index] <- score[index] + ibd.factor * (activ.posi - prev) / chrom.length
      }

      prev <- activ.posi

      activ[activ.min] <- min(activ[activ.min] +1, length(chr[[activ.min]]))
      activ.recom[activ.min] <- chr[[activ.min]][activ[activ.min]]
      activ.ursprung[activ.min] <-  origin[[activ.min]][activ[activ.min]]

    }
  }
  if(ibd.obs>0){
    ibd <- mean(score[1:ibd.obs])
  } else{
    ibd <- 0
  }
  if(hbd.obs>0){
    hbd <- mean(score[(ibd.obs+1):(ibd.obs+hbd.obs)])
  } else{
    hbd <- 0
  }

  return(c(ibd,hbd))
}





#' Approximate empirical kinship
#'
#' Function to compute empirical kinship for a set of individuals (not all pairs of individuals are evaluated)
#' @param population Population list
#' @param database1 First Groups of individuals to consider for the export
#' @param gen1 Quick-insert for database1 (vector of all generations to export)
#' @param cohorts1 Quick-insert for database1 (vector of names of cohorts to export)
#' @param database2 Second Groups of individuals to consider for the export
#' @param gen2 Quick-insert for database2 (vector of all generations to export)
#' @param cohorts2 Quick-insert for database2 (vector of names of cohorts to export)
#' @param ibd.obs Number of Individual pairs to sample for IBD estimation
#' @examples
#' data(ex_pop)
#' kinship.emp.fast(population=ex_pop,gen=2)
#' @return Empirical kinship matrix (IBD-based since Founders) per gen/database/cohort
#' @export
#'
kinship.emp.fast.between <- function(population=NULL, gen1=NULL, database1=NULL, cohorts1=NULL, gen2=NULL, database2=NULL, cohorts2=NULL,ibd.obs=50){
  if(length(gen1)>0 || length(database1)>0 || length(cohorts1)>0){
    database1 <- get.database(population, gen=gen1, database=database1, cohorts=cohorts1)
    animals1 <- list()
    for(index in 1:nrow(database1)){
      if(diff(database1[index,3:4])>=0){
        for(index2 in database1[index,3]:database1[index,4]){
          animals1[[length(animals1)+1]] <- population$breeding[[database1[index,1]]][[database1[[index,2]]]][[index2]]
        }
      }
    }
  }
  n1 <- length(animals1)

  if(length(gen2)>0 || length(database2)>0 || length(cohorts2)>0){
    database2 <- get.database(population, gen=gen2, database=database2, cohorts=cohorts2)
    animals2 <- list()
    for(index in 1:nrow(database2)){
      if(diff(database2[index,3:4])>=0){
        for(index2 in database2[index,3]:database2[index,4]){
          animals2[[length(animals2)+1]] <- population$breeding[[database2[index,1]]][[database2[[index,2]]]][[index2]]
        }
      }
    }
  }
  n2 <- length(animals2)

  if(n1==0 || n2==0){
    return("Input in both groups of animals required!")
  }

  chrom.length <- max(animals1[[1]][[1]])

  if(chrom.length=="removed"){
    return("Recombination history has been deleted!")
  }

  if(n1*n2 <= (ibd.obs)){
    ibd.obs <- n1*n2
    i1 <- rep(1:n1, length.out=ibd.obs)
    j1 <- sort(rep(1:n2, length.out=ibd.obs))
  } else{
    i1 <- sample(1:n1, ibd.obs, replace=if((ibd.obs)>n1){TRUE} else{FALSE})
    j1 <- sample(1:n2, ibd.obs, replace=if((ibd.obs)>n2){TRUE} else{FALSE})
  }


  score <- numeric(length(i1))
  for(index in 1:length(i1)){
    i <- i1[index]
    j <- j1[index]
    chr <- list()
    chr[[1]] <- animals1[[i]][[1]][-1]
    chr[[2]] <- animals1[[i]][[2]][-1]
    chr[[3]] <- animals2[[j]][[1]][-1]
    chr[[4]] <- animals2[[j]][[2]][-1]
    origin <- list()
    origin[[1]] <- animals1[[i]][[5]]
    origin[[2]] <- animals1[[i]][[6]]
    origin[[3]] <- animals2[[j]][[5]]
    origin[[4]] <- animals2[[j]][[6]]

    activ <- c(1,1,1,1)
    prev <- 0
    activ.recom <- c(chr[[1]][activ[1]], chr[[2]][activ[2]], chr[[3]][activ[3]], chr[[4]][activ[4]])
    activ.ursprung <- c(origin[[1]][activ[1]],origin[[2]][activ[2]],origin[[3]][activ[3]],origin[[4]][activ[4]])



    for(steps in 1:(length(c(chr[[1]], chr[[2]], chr[[3]], chr[[4]]))-3)){
      activ.min <- which.min(activ.recom)[1]
      activ.posi <- chr[[activ.min]][activ[activ.min]]

      if((activ.posi - prev)>0){
        #ibd <- length(unique(activ.ursprung)) # Nur vergleich des Neuen mit bisherigen Rechenzeiteffizienter!
        #ibd <- sum(!duplicated(activ.ursprung)) # Nur vergleich des Neuen mit bisherigen Rechenzeiteffizienter!
        #ibd.factor <- 1-ibd*0.25 + 0.25*(ibd==1)

        ibd.factor <- (sum(activ.ursprung[1] == activ.ursprung[3:4]) + sum(activ.ursprung[2] == activ.ursprung[3:4]))/4
        score[index] <- score[index] + ibd.factor * (activ.posi - prev) / chrom.length
      }

      prev <- activ.posi

      activ[activ.min] <- min(activ[activ.min] +1, length(chr[[activ.min]]))
      activ.recom[activ.min] <- chr[[activ.min]][activ[activ.min]]
      activ.ursprung[activ.min] <-  origin[[activ.min]][activ[activ.min]]

    }
  }
  if(ibd.obs>0){
    ibd <- mean(score[1:ibd.obs])
  } else{
    ibd <- 0
  }


  return(c(ibd))
}
