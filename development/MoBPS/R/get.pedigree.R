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

#' Derive pedigree
#'
#' Derive pedigree for selected individuals
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param founder.zero Parents of founders are displayed as "0" (default: TRUE)
#' @param raw Set to TRUE to not convert numbers into Sex etc.
#' @param id Set to TRUE to extract individual IDs
#' @examples
#' data(ex_pop)
#' get.pedigree(ex_pop, gen=2)
#' @return Pedigree-file for in gen/database/cohorts selected individuals
#' @export


get.pedigree <- function(population, database=NULL, gen=NULL, cohorts=NULL, founder.zero=TRUE,
                         raw=FALSE, id=FALSE, use.first.copy = FALSE){

  database <- get.database(population, gen, database, cohorts)

  n.animals <- sum(database[,4] - database[,3] +1)

  pedigree <- matrix(0, nrow=n.animals, ncol=3 + 6 * raw)
  rindex <- 1

  if(raw){
    for(row in 1:nrow(database)){
      animals <- database[row,]
      for(index in database[row,3]:database[row,4]){
        father <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[7]][1:3]
        mother <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[8]][1:3]

        if(use.first.copy){
          father = population$breeding[[father[1]]][[father[[2]]]][[father[[3]]]][[21]][1,]
          mother = population$breeding[[mother[1]]][[mother[[2]]]][[mother[[3]]]][[21]][1,]
        }
        pedigree[rindex,] <- c(database[row,1:2], index, father, mother)
        rindex <- rindex + 1
      }
    }
  } else{
    for(row in 1:nrow(database)){
      animals <- database[row,]
      for(index in database[row,3]:database[row,4]){
        father <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[7]]
        mother <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[8]]
        if(length(population$breeding[[father[1]]])>1){
          father_t <- population$breeding[[father[1]]][[father[2]+14]][[father[3]]]
        } else{
          father_t <- 0
        }
        if(length(population$breeding[[mother[1]]])>1){
          mother_t <- population$breeding[[mother[1]]][[mother[2]+14]][[mother[3]]]
        } else{
          mother_t <- 0
        }


        child_t <- population$breeding[[database[row,1]]][[database[row,2]+14]][[index]]
        pedigree[rindex,] <- c(child_t, father_t, mother_t)
        rindex <- rindex + 1
      }
    }
  }

  pedigree_id = pedigree

  if(!raw && !id){
    rindex <- 1
    for(row in 1:nrow(database)){
      animals <- database[row,]
      for(index in database[row,3]:database[row,4]){
        father <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[7]]
        mother <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[8]]

        if(use.first.copy){
          father = population$breeding[[father[1]]][[father[[2]]]][[father[[3]]]][[21]][1,]
          mother = population$breeding[[mother[1]]][[mother[[2]]]][[mother[[3]]]][[21]][1,]
        }

        father_t <- paste(if(father[2]==1) "M" else "F", father[3], "_", father[1], sep="")
        mother_t <- paste(if(mother[2]==1) "M" else "F", mother[3], "_", mother[1], sep="")
        child_t <- paste(if(database[row,2]==1) "M" else "F", index, "_", database[row,1], sep="")
        pedigree[rindex,] <- c(child_t, father_t, mother_t)
        rindex <- rindex + 1
      }
    }
  }



  if(!raw){
    colnames(pedigree) <- c("offspring", "father", "mother")
  } else{
    colnames(pedigree) <- c("offspring.gen","offspring.sex","offspring.nr",
                            "father.gen", "father.sex","father.nr",
                            "mother.gen","mother.sex","mother.nr")
  }


  if(founder.zero && !raw){
    set0 <- which(pedigree_id[,1]==pedigree_id[,2])
    if(length(set0)>0){
      pedigree[set0,2] <- "0"
    }
    set0 <- which(pedigree_id[,1]==pedigree_id[,3])
    if(length(set0)>0){
      pedigree[set0,3] <- "0"
    }
  }
  return(pedigree)
}

#' Derive pedigree including grandparents
#'
#' Derive pedigree for selected individuals including grandparents
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param shares Determine actual inherited shares of grandparents
#' @param founder.zero Parents of founders are displayed as "0" (default: TRUE)
#' @param raw Set to TRUE to not convert numbers into Sex etc.
#' @examples
#' data(ex_pop)
#' get.pedigree2(ex_pop, gen=3)
#' @return Pedigree-file (grandparents) for in gen/database/cohorts selected individuals
#' @export

get.pedigree2 <- function(population, database=NULL, gen=NULL, cohorts=NULL, shares=FALSE, founder.zero=TRUE, raw=FALSE){

  database <- get.database(population, gen, database, cohorts)

  n.animals <- sum(database[,4] - database[,3] +1)


  if(shares){
    cols <- 9 + 10*raw
  } else{
    cols <- 5 + 10*raw
  }
  pedigree <- matrix(0, nrow=n.animals, ncol=cols)
  rindex <- 1

  if(raw){
    for(row in 1:nrow(database)){
      animals <- database[row,]
      for(index in database[row,3]:database[row,4]){
        father <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[7]][1:3]
        mother <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[8]][1:3]
        grandf1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[7]][1:3]
        grandf1share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandm1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[8]][1:3]
        grandm1share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandf2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[7]][1:3]
        grandf2share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]
        grandm2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[8]][1:3]
        grandm2share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]

        grandf1_t <- paste(if(grandf1[2]==1) "M" else "F", grandf1[3], "_", grandf1[1], sep="")
        grandm1_t <- paste(if(grandm1[2]==1) "M" else "F", grandm1[3], "_", grandm1[1], sep="")
        grandf2_t <- paste(if(grandf2[2]==1) "M" else "F", grandf2[3], "_", grandf2[1], sep="")
        grandm2_t <- paste(if(grandm2[2]==1) "M" else "F", grandm2[3], "_", grandm2[1], sep="")
        child_t <- paste(if(database[row,2]==1) "M" else "F", index, "_", gen, sep="")

        father_t <- paste(if(father[2]==1) "M" else "F", father[3], "_", father[1], sep="")
        mother_t <- paste(if(mother[2]==1) "M" else "F", mother[3], "_", mother[1], sep="")
        child_t <- paste(if(database[row,2]==1) "M" else "F", index, "_", database[row,1], sep="")
        if(shares){

          if(length(grandf1share)==0){
            grandf1share <- NA
          }
          if(length(grandm1share)==0){
            grandm1share <- NA
          }
          if(length(grandf2share)==0){
            grandf2share <- NA
          }
          if(length(grandm2share)==0){
            grandm2share <- NA
          }

          pedigree[rindex,] <-c(animals[1:2], index,   grandf1, grandm1, grandf2, grandm2, grandf1share, grandm1share, grandf2share, grandm2share)

          colnames(pedigree) <- c("offspring.gen","offspring.sex","offspring.nr",
                                  "paternal grandfather.gen","paternal grandfather.sex","paternal grandfather.nr",
                                  "paternal grandmother.gen", "paternal grandmother.sex", "paternal grandmothe.nr",
                                  "maternal grandfather.gen",  "maternal grandfather.sex", "maternal grandfather.nr",
                                  "maternal grandmother.gen", "maternal grandmother.sex", "maternal grandmother.nr",
                                  "paternal grandfather.share", "paternal grandmother.share",
                                  "maternal grandfather.share", "maternal grandmother.share")
        } else{
          pedigree[rindex,] <- c(animals[1:2], index,   grandf1, grandm1, grandf2, grandm2)

          colnames(pedigree) <- c("offspring.gen","offspring.sex","offspring.nr",
                                  "paternal grandfather.gen","paternal grandfather.sex","paternal grandfather.nr",
                                  "paternal grandmother.gen", "paternal grandmother.sex", "paternal grandmothe.nr",
                                  "maternal grandfather.gen",  "maternal grandfather.sex", "maternal grandfather.nr",
                                  "maternal grandmother.gen", "maternal grandmother.sex", "maternal grandmother.nr")
        }
        rindex <- rindex + 1
      }
    }
  } else{
    for(row in 1:nrow(database)){
      animals <- database[row,]
      for(index in database[row,3]:database[row,4]){
        father <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[7]][1:3]
        mother <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[8]][1:3]
        grandf1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[7]][1:3]
        grandf1share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandm1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[8]][1:3]
        grandm1share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandf2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[7]][1:3]
        grandf2share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]
        grandm2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[8]][1:3]
        grandm2share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]

        grandf1_t <- paste(if(grandf1[2]==1) "M" else "F", grandf1[3], "_", grandf1[1], sep="")
        grandm1_t <- paste(if(grandm1[2]==1) "M" else "F", grandm1[3], "_", grandm1[1], sep="")
        grandf2_t <- paste(if(grandf2[2]==1) "M" else "F", grandf2[3], "_", grandf2[1], sep="")
        grandm2_t <- paste(if(grandm2[2]==1) "M" else "F", grandm2[3], "_", grandm2[1], sep="")
        child_t <- paste(if(database[row,2]==1) "M" else "F", index, "_", database[row,1], sep="")

        if(shares){
          if(length(grandf1share)==0){
            grandf1share <- NA
          }
          if(length(grandm1share)==0){
            grandm1share <- NA
          }
          if(length(grandf2share)==0){
            grandf2share <- NA
          }
          if(length(grandm2share)==0){
            grandm2share <- NA
          }

          pedigree[rindex,] <-  c(child_t, grandf1_t, grandm1_t, grandf2_t, grandm2_t, grandf1share, grandm1share, grandf2share, grandm2share)

          colnames(pedigree) <- c("offspring",
                                  "paternal grandfather",
                                  "paternal grandmother",
                                  "maternal grandfather",
                                  "maternal grandmother",
                                  "paternal grandfather.share", "paternal grandmother.share",
                                  "maternal grandfather.share", "maternal grandmother.share")
        } else{
          pedigree[rindex,] <- c(child_t, grandf1_t, grandm1_t, grandf2_t, grandm2_t)
          colnames(pedigree) <- c("offspring",
                                  "paternal grandfather",
                                  "paternal grandmother",
                                  "maternal grandfather",
                                  "maternal grandmother")
        }
        rindex <- rindex + 1
      }
    }
  }


  if(founder.zero){
    pedi1 <- get.pedigree(population, database = database)
  }

  if(founder.zero){
    set0 <- which(pedigree[,1]==pedigree[,2] | pedi1[,2]==pedigree[,2])
    if(length(set0)>0){
      pedigree[set0,2] <- "0"
    }
    set0 <- which(pedigree[,1]==pedigree[,3]| pedi1[,2]==pedigree[,3])
    if(length(set0)>0){
      pedigree[set0,3] <- "0"
    }
  }
  if(founder.zero){
    set0 <- which(pedigree[,1]==pedigree[,4] | pedi1[,3]==pedigree[,4])
    if(length(set0)>0){
      pedigree[set0,4] <- "0"
    }
    set0 <- which(pedigree[,1]==pedigree[,5] | pedi1[,3]==pedigree[,5])
    if(length(set0)>0){
      pedigree[set0,5] <- "0"
    }
  }
  return(pedigree)
}

#' Derive pedigree parents and grandparents
#'
#' Derive pedigree for selected individuals including parents/grandparents
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param founder.zero Parents of founders are displayed as "0" (default: TRUE)
#' @param raw Set to TRUE to not convert numbers into Sex etc.
#' @param id Set to TRUE to extract individual IDs
#' @examples
#' data(ex_pop)
#' get.pedigree3(ex_pop, gen=3)
#' @return Pedigree-file (parents + grandparents) for in gen/database/cohorts selected individuals
#' @export
#'
get.pedigree3 <- function(population, database=NULL, gen=NULL, cohorts=NULL, founder.zero=TRUE, id=FALSE, raw=FALSE){

  database <- get.database(population, gen, database, cohorts)

  n.animals <- sum(database[,4] - database[,3] +1)
  pedigree <- matrix(0, nrow=n.animals, ncol= 7 + 14*raw)
  rindex <- 1

  if(raw){
    for(row in 1:nrow(database)){
      animals <- database[row,]
      for(index in database[row,3]:database[row,4]){
        father <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[7]][1:3]
        mother <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[8]][1:3]
        grandf1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[7]][1:3]
        grandf1share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandm1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[8]][1:3]
        grandm1share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandf2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[7]][1:3]
        grandf2share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]
        grandm2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[8]][1:3]
        grandm2share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]


        pedigree[rindex,] <- c(database[row,1:2], index, father, mother, grandf1, grandm1, grandf2, grandm2)
        rindex <- rindex + 1
      }
    }
    colnames(pedigree) <- c("offspring.gen","offspring.sex","offspring.nr",
                            "father.gen", "father.sex","father.nr",
                            "mother.gen","mother.sex","mother.nr",
                            "paternal grandfather.gen","paternal grandfather.sex","paternal grandfather.nr",
                            "paternal grandmother.gen", "paternal grandmother.sex", "paternal grandmothe.nr",
                            "maternal grandfather.gen",  "maternal grandfather.sex", "maternal grandfather.nr",
                            "maternal grandmother.gen", "maternal grandmother.sex", "maternal grandmother.nr")
  } else {
    for(row in 1:nrow(database)){
      animals <- database[row,]
      for(index in database[row,3]:database[row,4]){
        father <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[7]]
        mother <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[8]]
        grandf1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[7]]
        grandf1share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandm1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[8]]
        grandm1share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandf2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[7]]
        grandf2share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]
        grandm2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[8]]
        grandm2share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]

        if(length(population$breeding[[father[1]]])>1){
          f1_t <- population$breeding[[father[1]]][[father[2]+14]][[father[3]]]
        } else{
          f1_t <- 0
        }

        if(length(population$breeding[[mother[1]]])>1){
          m1_t <- population$breeding[[mother[1]]][[mother[2]+14]][[mother[3]]]
        } else{
          m1_t <- 0
        }

        if(length(population$breeding[[grandf1[1]]])>1){
          grandf1_t <- population$breeding[[grandf1[1]]][[grandf1[2]+14]][[grandf1[3]]]
        } else{
          grandf1_t <- 0
        }

        if(length(population$breeding[[grandf2[1]]])>1){
          grandf2_t <- population$breeding[[grandf2[1]]][[grandf2[2]+14]][[grandf2[3]]]
        } else{
          grandf2_t <- 0
        }

        if(length(population$breeding[[grandm1[1]]])>1){
          grandm1_t <- population$breeding[[grandm1[1]]][[grandm1[2]+14]][[grandm1[3]]]
        } else{
          grandm1_t <- 0
        }

        if(length(population$breeding[[grandm2[1]]])>1){
          grandm2_t <- population$breeding[[grandm2[1]]][[grandm2[2]+14]][[grandm2[3]]]
        } else{
          grandm2_t <- 0
        }

        child_t <- population$breeding[[database[row,1]]][[database[row,2]+14]][[index]]

        pedigree[rindex,] <- c(child_t, f1_t, m1_t, grandf1_t, grandm1_t, grandf2_t, grandm2_t)
        rindex <- rindex + 1
      }
    }
    colnames(pedigree) <- c("offspring", "father", "mother", "paternal grandfather", "paternal grandmother", "maternal grandfather", "maternal grandmother")

  }

  pedigree_id = pedigree
  if(!raw && !id){
    rindex = 1
    for(row in 1:nrow(database)){
      animals <- database[row,]
      for(index in database[row,3]:database[row,4]){
        father <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[7]]
        mother <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[8]]
        grandf1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[7]]
        grandf1share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandm1 <- population$breeding[[father[[1]]]][[father[[2]]]][[father[[3]]]][[8]]
        grandm1share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[19]]
        grandf2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[7]]
        grandf2share <- population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]
        grandm2 <- population$breeding[[mother[[1]]]][[mother[[2]]]][[mother[[3]]]][[8]]
        grandm2share <- 1 - population$breeding[[database[row,1]]][[database[row,2]]][[index]][[20]]

        f1_t <- paste(if(father[2]==1) "M" else "F", father[3], "_", father[1], sep="")
        m1_t <- paste(if(mother[2]==1) "M" else "F", mother[3], "_", mother[1], sep="")
        grandf1_t <- paste(if(grandf1[2]==1) "M" else "F", grandf1[3], "_", grandf1[1], sep="")
        grandm1_t <- paste(if(grandm1[2]==1) "M" else "F", grandm1[3], "_", grandm1[1], sep="")
        grandf2_t <- paste(if(grandf2[2]==1) "M" else "F", grandf2[3], "_", grandf2[1], sep="")
        grandm2_t <- paste(if(grandm2[2]==1) "M" else "F", grandm2[3], "_", grandm2[1], sep="")
        child_t <- paste(if(database[row,2]==1) "M" else "F", index, "_", database[row,1], sep="")

        father_t <- paste(if(father[2]==1) "M" else "F", father[3], "_", father[1], sep="")
        mother_t <- paste(if(mother[2]==1) "M" else "F", mother[3], "_", mother[1], sep="")
        child_t <- paste(if(database[row,2]==1) "M" else "F", index, "_", database[row,1], sep="")

        pedigree[rindex,] <- c(child_t, f1_t, m1_t, grandf1_t, grandm1_t, grandf2_t, grandm2_t)
        rindex <- rindex + 1
      }
    }
    colnames(pedigree) <- c("offspring", "father", "mother", "paternal grandfather", "paternal grandmother", "maternal grandfather", "maternal grandmother")
  }

  if(!raw){

    if(founder.zero){
      set0 <- which(pedigree_id[,2]==pedigree_id[,4])
      if(length(set0)>0){
        pedigree[set0,4] <- 0
      }
      set0 <- which(pedigree_id[,2]==pedigree_id[,5])
      if(length(set0)>0){
        pedigree[set0,5] <- 0
      }
    }
    if(founder.zero){
      set0 <- which(pedigree_id[,3]==pedigree_id[,6])
      if(length(set0)>0){
        pedigree[set0,6] <- 0
      }
      set0 <- which(pedigree_id[,3]==pedigree_id[,7])
      if(length(set0)>0){
        pedigree[set0,7] <- 0
      }
    }

    if(founder.zero){
      set0 <- which(pedigree_id[,1]==pedigree_id[,2])
      if(length(set0)>0){
        pedigree[set0,2] <- 0
      }
      set0 <- which(pedigree_id[,1]==pedigree_id[,3])
      if(length(set0)>0){
        pedigree[set0,3] <- 0
      }
    }

  }



  return(pedigree)
}
