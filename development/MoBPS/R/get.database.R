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

#' gen/database/cohorts conversion
#'
#' Function to derive a database based on gen/database/cohorts
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param avoid.merging Set to TRUE to avoid different cohorts to be merged in a joint group when possible
#' @param id Individual IDs to search/collect in the database
#' @param id.all.copy Set to TRUE to show all copies of an individual in the database (default: FALSE)
#' @param id.last Set to TRUE to use the last copy of an individual for the database (default: FALSE - pick first copy)
#' @examples
#' data(ex_pop)
#' get.database(ex_pop, gen=2)
#' @return Combine gen/database/cohorts to a joined database
#' @export


get.database<- function(population, gen=NULL, database=NULL, cohorts=NULL, avoid.merging=FALSE, id=NULL, id.all.copy=FALSE, id.last=FALSE){

  if(length(id)>0 && (length(gen)>0 || length(database)>0 || length(cohorts) > 0 )){
    stop("You can either prove IDs or gen/database/cohorts")
  }

  if(length(id)>0){
    if(ncol(population$info$cohorts)<10){

      min_max <- NULL
      coh <- get.cohorts(population)
      for(index in coh){
        temp1 <- get.id(population, cohorts = index)
        min_max <- rbind(min_max, c(min(temp1), max(temp1)))

      }
    }

    database <- matrix(0, nrow=length(id), ncol=4)
    k <- 1
    for(index in id){
      poss <- which((index <= min_max[,2]) & (index >= min_max[,1]))

      ids <- get.id(population, cohorts = coh[poss])
      db <- get.database(population, cohorts=coh[poss])

      if(id.last){
        temp1 <- which(ids==index)
        activ <- temp1[length(temp1)]
      } else{
        activ <- which(ids==index)[1]
      }


      n.animals <- c(0,cumsum(db[,4]-db[,3]+1))
      which_coh <- which(activ > n.animals)
      which_coh <- max(which_coh)
      which_nr <- activ - n.animals[which_coh] + db[which_coh,3] - 1
      database[k,] <- c(db[which_coh,1:2], which_nr, which_nr)
      k <- k + 1
    }



    if(!id.all.copy){
      database <- get.database(population, database=database)
      return(database)
    } else{
      database_list <- list()

      for(index in 1:nrow(database)){
        database_list[[index]] <- population$breeding[[database[index,1]]][[database[index,2]]][[database[index,3]]][[21]]
      }

      return(database_list)

    }

  }


  if(length(gen)>0){
    database_gen <- cbind(rep(gen,each=2), rep(1:2, length(gen)))
    if(length(database)>0 && ncol(database)==4){
      database_gen <- cbind(database_gen,1,population$info$size[database_gen])
    }
    database <- rbind(database, database_gen)
  }
  if(length(database)>0 && !is.matrix(database)){
    database <- matrix(database, nrow=1)
  }
  if(length(database)>0 && ncol(database)==2){
    start <- end <- numeric(nrow(database))
    for(index in 1:nrow(database)){
      start[index] <- 1
      end[index] <- population$info$size[database[index,1], database[index,2]]
    }
    database <- cbind(database, start, end)
  }
  if(length(database)>0 && ncol(database)==3){
    database <- cbind(database, database[,3])
  }



  if(length(cohorts)>0){
    database2 <- matrix(NA, nrow=length(cohorts)*3, ncol=4)
    added = 0
    remove = NULL
    ncoh = length(cohorts)
    for(index in 1:(length(cohorts)*3)){
      if(length(cohorts)< index){
        break
      }
      row <- which(population$info$cohorts[,1]==cohorts[index])[1]
      if(is.na(row)){


        candidates = paste0(cohorts[index], c("_M", "_F"))
        if(sum(population$info$cohorts[,1]==candidates[1])>0 && sum(population$info$cohorts[,1]==candidates[2])>0){
          cohorts = c(cohorts, candidates)
          added = added + 2
          remove = c(remove, index)
          warning(paste0("Cohort ", cohorts[index], " is not available. Instead use ",  cohorts[index], "_M & ",  cohorts[index], "_F!"))
        } else{
          warning(paste0("Cohort ", cohorts[index], " is not available when constructing the individuals database!"))
        }

      }
      gen <- as.numeric(population$info$cohorts[row,2])
      sex <- 1 + (as.numeric(population$info$cohorts[row,4])>0)
      first <- as.numeric(population$info$cohorts[row,5 + sex])
      last <- first + as.numeric(population$info$cohorts[row,2 + sex]) - 1
      database2[index,] <- c(gen,sex,first,last)
    }
    database2 = database2[1:(ncoh + added),,drop=FALSE]
    if(length(remove)>0){
      database2 = database2[-remove,, drop=FALSE]
    }

    if(sum(is.na(database2))>0){
      warning("Cohort-name is not available! \nCheck cohort names (in particular for added '_F' and '_M') / get.cohorts()!")
      database2 <- database2[!is.na(database2[,1]), ,drop=FALSE]
    }
    database <- rbind(database, database2)
  }

  keep <- database[,3]<=database[,4]
  database <- database[keep,,drop=FALSE]

  if(length(database)>0 && nrow(database)>1){
    order <- sort(database[,1]*1e12 + database[,2]*1e7 + database[,3], index.return=TRUE)$ix
    database <- database[order,,drop=FALSE]
    first_same <- 1
    first_index <- 1
    not_first = FALSE
    for(index in 2:nrow(database)){
      if(database[first_index,1]!=database[index,1] || database[first_index,2]!=database[index,2]){
        first_index <- which(database[index,1]==database[,1])[1]
        first_same <- database[first_index,1]
        not_first = FALSE
      } else{
        if(database[(index-1),1]!=0){
          if(database[index-1,1]==database[index,1] & database[index-1,2] == database[index,2]){
            checks <- (index-1)
          } else{
            checks <- NULL
          }
        } else{
          if(not_first){
            checks <- (which(database[first_index:(index-1),1]==database[index,1] & database[first_index:(index-1),2] == database[index,2])) + first_index - 1
          } else{
            checks <- first_index

          }
        }
        if(!avoid.merging){
          for(index2 in checks){
            if(database[index,3] <= (database[index2,4]+1)){
              database[index2,4] <- max(database[index2,4], database[index,4])
              database[index,] <- 0
            } else{
              not_first = TRUE
            }

          }
        }

      }

    }
    database <- database[database[,1]!=0,,drop=FALSE]
  }

  database <- unique(database)

  return(database)
}
