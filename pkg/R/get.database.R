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

#' gen/database/cohorts conversion
#'
#' Function to derive a database based on gen/database/cohorts
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param avoid.merging Set to TRUE to avoid different cohorts to be merged in a joint group when possible
#' @param id Individual IDs to search/collect in the database
#' @param db.names MoPBS internal names (SexNr_Generation)
#' @param keep.order To not change order of individuals when ids are provided (default: FALSE)
#' @param id.all.copy Set to TRUE to show all copies of an individual in the database (default: FALSE)
#' @param id.last Set to TRUE to use the last copy of an individual for the database (default: FALSE - pick first copy)
#' @param class Only include individuals of the following classes in the database (can also be vector with multiple classes; default: ALL)
#' @param per.individual Set TRUE to obtain a database with one row per individual instead of concatenating (default: FALSE)
#' @param sex.filter Set to 1 to only include males and set 2 to only include females in database (default: 0)
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' data(ex_pop)
#' get.database(ex_pop, gen=2)
#' @return Matrix with combined gen/database/cohorts
#' @export


get.database<- function(population, gen=NULL, database=NULL, cohorts=NULL, avoid.merging=FALSE,
                        per.individual = FALSE, id=NULL, db.names = NULL, id.all.copy=FALSE, id.last=FALSE,
                        keep.order = FALSE, class = NULL, verbose = TRUE,
                        sex.filter = 0){

  if(length(id)>0 && (length(gen)>0 || length(database)>0 || length(cohorts) > 0 )){
    stop("You can either prove IDs or gen/database/cohorts")
  }

  if(length(db.names)>0){

    sex = substr(db.names, start = 1, stop = 1)
    rest = substr(db.names, start = 2, stop = 100)
    rest_split = strsplit(rest, split = "_")

    database_names = cbind(0, as.numeric(sex == "F") + 1, 0,0)
    for(index in 1:length(rest_split)){

      database_names[index,c(1,3,4)] = rest_split[[index]][c(2,1,1)]

    }

    storage.mode(database_names) = "numeric"

    database = rbind(get.database(population, database = database), database_names)

    database = get.database(population, gen=gen, database=database, cohorts=cohorts, avoid.merging=avoid.merging,
                            per.individual = per.individual, id=id, db.names = NULL, id.all.copy=id.all.copy, id.last=id.last,
                            keep.order = keep.order, class = class, verbose = verbose)
    return(database)

  }

  if(length(id)>0){
    if(is.character(id[1])){
      id = as.numeric(id)
    }

    coh <- get.cohorts(population)

    if(ncol(population$info$cohorts)<10){
      min_max <- NULL
      for(index in coh){
        temp1 <- get.id(population, cohorts = index)
        min_max <- rbind(min_max, c(min(temp1), max(temp1)))

      }
    } else{
      min_max = (population$info$cohorts[,10:11,drop=FALSE])
      storage.mode(min_max) = "numeric"
    }

    database <- matrix(0, nrow=length(id), ncol=4)
    k <- 1

    poss_prior = NULL
    for(index in id){
      poss <- which((index <= min_max[,2]) & (index >= min_max[,1]))

      if(!(length(poss)==length(poss_prior) && prod(poss==poss_prior)==1)){
        ids <- get.id(population, cohorts = coh[poss])
        db <- get.database(population, cohorts=coh[poss])
        poss_prior = poss
      }



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

    if(keep.order){
      if(per.individual){
        database <- get.database(population, database=database, per.individual = per.individual)
      }
      return(database)
    }

    if(!id.all.copy){
      database <- get.database(population, database=database, per.individual = per.individual)
      return(database)
    } else{
      database_list <- list()

      for(index in 1:nrow(database)){
        database_list[[index]] <- population$breeding[[database[index,1]]][[database[index,2]]][[database[index,3]]][[21]]
      }

      return(database_list)

    }

  }

  if(length(class)>0 && length(gen)==0 && length(database)==0 && length(cohorts)==0){
    gen = 1:get.ngen(population)
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
          # if(verbose) cat(paste0("Cohort ", cohorts[index], " is not available. Instead use ",  cohorts[index], "_M & ",  cohorts[index], "_F!\n"))
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

   # database = rbind(pedigree.database[300:350,], cbind(5,2,1,400))


    tmp  = sort(database[,1]*1e12 + database[,2]*1e7 + database[,3], index.return=TRUE)

    order <- tmp$ix
    database <- database[order,,drop=FALSE]

    if(nrow(database)>1){
      tmp2 = ((database[-nrow(database),4] + 1) == database[-1,3]) & ((database[-nrow(database),2]) == database[-1,2]) & ((database[-nrow(database),1]) == database[-1,1])
    } else{
      tmp2 = rep(FALSE, length(tmp)-1)
    }

    # first very basic merging
    first = which(c(TRUE, diff(tmp$x)!=1 | !tmp2 ))
    last = which(c(diff(tmp$x)!=1 | !tmp2 , TRUE))
    database = cbind(database[first,1:3,drop = FALSE], database[last,4,drop = FALSE])

    if(nrow(database) > 1){
      first_same <- 1
      first_index <- 1
      not_first = FALSE
      for(index in 2:nrow(database)){
        if(database[first_index,1]!=database[index,1] || database[first_index,2]!=database[index,2]){
          first_index <- which(database[index,1]==database[,1] & database[index,2]==database[,2])[1]
          first_same <- database[first_index,1]
          not_first = FALSE
        } else{
          if(database[(index-1),1]!=0){
            if(database[index-1,1]==database[index,1] & database[index-1,2] == database[index,2]){
              checks <- (index-1)
            } else{
              checks <- NULL
            }
          } else if(database[(index-2),1]!=0 && first_index <= index-2){
            if(database[(index-2),1]!=0){
              if(database[index-2,1]==database[index,1] & database[index-2,2] == database[index,2]){
                checks <- (index-2)
              } else{
                checks <- NULL
              }
            }
          } else if(database[(index-3),1]!=0 && first_index <= index-3){
            if(database[(index-3),1]!=0){
              if(database[index-3,1]==database[index,1] & database[index-3,2] == database[index,2]){
                checks <- (index-3)
              } else{
                checks <- NULL
              }
            }
          } else{
            if(not_first){
              checks <- (which(database[first_index:(index-1),1]==database[index,1] & database[first_index:(index-1),2] == database[index,2])) + first_index - 1
            } else{
              checks <- first_index

            }
          }
          if(!avoid.merging){

            checks = checks[database[index,3] <= (database[checks,4] + 1)]

            if(length(checks)>0){
              database[checks,4] <- max(database[checks,4], database[index,4])
              database[index,] <- 0
            }
            not_first = TRUE
          }

        }

      }
    }

    database <- database[database[,1]!=0,,drop=FALSE]
  }

  if(sex.filter != 0){
    database = database[database[,2] == sex.filter]
  }

  database <- unique(database)

  if(length(class)>0){

    class_db = get.class(population, database = database)

    temp_db = matrix(0, nrow = length(class_db), ncol = 3)
    current = 1
    for(index in 1:nrow(database)){
      n_indi = database[index,4] - database[index,3] + 1
      temp_db[current:(current + n_indi -1),1] = database[index,1]
      temp_db[current:(current + n_indi -1),2] = database[index,2]
      temp_db[current:(current + n_indi -1),3] = database[index,3]:database[index,4]
      current = current + n_indi
    }

    keep = rep(FALSE, length(class_db))
    for(index in 1:length(class)){
      keep[class_db == class[index]] = TRUE
    }

    temp_db = temp_db[keep,,drop=FALSE]

    database = get.database(population, database = temp_db)
  }

  if(per.individual){
    n_animal = sum(database[,4] - database[,3] + 1)

    database_full = matrix(0, nrow = n_animal, ncol = 4)

    so_far = 0
    for(index in 1:nrow(database)){
      consider = database[index,3]:database[index,4]
      database_full[1:length(consider) + so_far,] = cbind(database[index,1], database[index,2],consider , consider )
      so_far = so_far + length(consider)
    }

    database = database_full
  }

  return(database)
}
