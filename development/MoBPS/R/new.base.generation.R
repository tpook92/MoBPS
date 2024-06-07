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

#' Set new base generation
#'
#' Function to set a new base generation for the population
#' @param population Population list
#' @param base.gen Vector containing all new base generations
#' @param base.database Matrix containing all database entries to be used as new base generation
#' @param base.cohorts Vector containing all cohorts to be used as new base generations
#' @param delete.previous.gen Delete all data before base.gen (default: FALSE)
#' @param delete.breeding.totals Delete all breeding totals before base.gen (default: FALSE)
#' @param delete.bve.data Deleta all previous bve data (default: FALSE)
#' @param add.chromosome.ends Add chromosome ends as recombination points
#' @param founder.pool AAA
#' @examples
#' data(ex_pop)
#' ex_pop <- new.base.generation(ex_pop, base.gen=2)
#' @return Population-List with mutated marker for the selected individual
#' @export

new.base.generation <- function(population, base.gen=NULL,
                                base.database = NULL,
                                base.cohorts = NULL,
                                delete.previous.gen=FALSE, delete.breeding.totals=FALSE,
                                delete.bve.data=FALSE, add.chromosome.ends=TRUE,
                                founder.pool = 1){

  if (requireNamespace("miraculix", quietly = TRUE)) {
    codeOriginsU <- miraculix::codeOrigins
    decodeOriginsU <- miraculix::decodeOrigins
  } else{
    codeOriginsU <- codeOriginsR
    decodeOriginsU <- decodeOriginsR
  }
  if(length(population$info$miraculix)>0 && population$info$miraculix){
    miraculix <- TRUE
  } else{
    miraculix <- FALSE
  }

  base.database.full = get.database(population, gen = base.gen, database = base.database, cohorts = base.cohorts)

  base.gen = unique(base.database.full[,1])

  if(length(base.gen)==0){
    base.gen <- length(population$breeding)
  }

  for(gen in base.gen){
    take <- which(population$info$origin.gen==gen)
    if(length(take)==1){
      origin_code <- population$info$origin.gen[take]
    } else{
      if(population$info$miraculix){
        if(length(population$info$origin.gen)<64){
          population$info$origin.gen <- c(population$info$origin.gen, as.integer(gen))
          origin_code <- length(population$info$origin.gen)
        } else{
          warning("To many origin generation!")
          warning("Delete second lowest origin.gen")
          switch_gen <- sort(population$info$origin.gen, index.return=TRUE)$ix[2]
          population$info$origin.gen[switch_gen] <- as.integer(gen)
          origin_code <- switch_gen
        }
      } else{
        if(length(population$info$origin.gen)<32){
          population$info$origin.gen <- c(population$info$origin.gen, as.integer(gen))
          origin_code <- length(population$info$origin.gen)
        } else{
          warning("To many origin generation!")
          warning("Delete second lowest origin.gen")
          switch_gen <- sort(population$info$origin.gen, index.return=TRUE)$ix[2]
          population$info$origin.gen[switch_gen] <- as.integer(gen)
          origin_code <- switch_gen
        }
      }
    }



    if(length(base.database) > 0 || length(base.cohorts)>0){
      skip = TRUE

    } else{
      skip = FALSE

    }

    for(sex in 1:2){
      if(length(population$breeding[[gen]][[sex]])>0){
        for(nr in 1:length(population$breeding[[gen]][[sex]])){

          if(skip){

            checks = (base.database.full[,1] == gen) + (base.database.full[,2] == sex) + (base.database.full[,3] <= nr) + (base.database.full[,4] >= nr)

            if(max(checks)<4){
              next
            }
          }

          if(miraculix){
            population$breeding[[gen]][[sex]][[nr]][[9]] <- miraculix::computeSNPS(population, gen, sex, nr, what="haplo", output_compressed=TRUE)
          } else{
            snps <- compute.snps(population, gen, sex, nr, decodeOriginsU=decodeOriginsU)
            population$breeding[[gen]][[sex]][[nr]][[9]] <- snps[1,]
            population$breeding[[gen]][[sex]][[nr]][[10]] <- snps[2,]
          }

          population$breeding[[gen]][[sex]][[nr]][[1]] <- c(0, sum(population$info$length))
          population$breeding[[gen]][[sex]][[nr]][[2]] <- c(0, sum(population$info$length))
          if(add.chromosome.ends==TRUE){
            population$breeding[[gen]][[sex]][[nr]][[1]] <- population$info$length.total
            population$breeding[[gen]][[sex]][[nr]][[2]] <- population$info$length.total
          }
          population$breeding[[gen]][[sex]][[nr]][[3]] <- numeric(0)
          population$breeding[[gen]][[sex]][[nr]][[4]] <- numeric(0)
          population$breeding[[gen]][[sex]][[nr]][[5]] <- codeOriginsU(matrix(c(origin_code, sex, nr, 1),nrow=(length(population$breeding[[gen]][[sex]][[nr]][[1]])-1), ncol=4, byrow=TRUE))

          population$breeding[[gen]][[sex]][[nr]][[6]] <- codeOriginsU(matrix(c(origin_code, sex, nr, 2),nrow=(length(population$breeding[[gen]][[sex]][[nr]][[2]])-1), ncol=4, byrow=TRUE))

          population$breeding[[gen]][[36+sex]][nr] <- founder.pool
        }

      }

    }
  }

  population$info$founder_pools = unique(c(population$info$founder_pools, founder.pool))
  population$info$founder_multi = if(length(population$info$founder_pools)>1){TRUE} else{FALSE}

  if(delete.previous.gen){
    for(index in 1:(min(base.gen)-1)){
      population$breeding[[index]] <- "deleted"
    }
  }
  if(delete.breeding.totals){
    population$info$breeding.totals <- list()
  }
  if(delete.bve.data){
    population$info$bve.data <- list()
  }
  return(population)
}
