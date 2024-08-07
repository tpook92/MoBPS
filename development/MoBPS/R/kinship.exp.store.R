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

#' Derive expected kinship
#'
#' Function to derive expected kinship
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param depth.pedigree Depth of the pedigree in generations
#' @param start.kinship Relationship matrix of the individuals in the first considered generation
#' @param elements Vector of individuals from the database to include in pedigree matrix
#' @param mult Multiplicator of kinship matrix (default: 1; set to 2 for a pedigree relationship matrix)
#' @param storage.save Lower numbers will lead to less memory but slightly higher computing time (default: 1.5, min: 1)
#' @param elements.copy Set to TRUE to automatically remove duplicated entries for individuals/database from the output matrix (default: FALSE)
#' @param include.error Set to FALSE to ignore/correct any errors in the pedigree
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' data(ex_pop)
#' kinship <- kinship.exp(population=ex_pop, gen=2)
#' @return Pedigree-based kinship matrix for in gen/database/cohort selected individuals
#' @export
#'

kinship.exp <- function(population, gen=NULL, database=NULL, cohorts=NULL, depth.pedigree=7,
                               start.kinship=NULL,
                               elements = NULL,
                               mult = 1,
                               storage.save=1.5,
                               verbose=TRUE,
                        include.error = TRUE,
                        elements.copy = FALSE){

  t1 <- as.numeric(Sys.time())
  int_mult <- as.integer(2^29)
  int_mult2 <- as.integer(2^28)
  database <- get.database(population, gen=gen, database=database, cohorts=cohorts)

  n.animals <- sum(diff(t(database[,3:4, drop=FALSE]))+1)

  if(length(elements)>0){

    if(max(elements)> n.animals){
      stop("kinship emp number of individuals does not match!")
    }

    cumorder <- cumsum(c(1,diff(t(database[,3:4, drop=FALSE]))+1))
    activ_database <- rep(FALSE, nrow(database))
    elements_new <- elements
    for(index in 1:nrow(database)){
      if(sum(intersect(elements, cumorder[index]:(cumorder[index+1]-1)))>0){
        activ_database[index] <- TRUE
      } else{
        elements_new[elements>cumorder[index]] <- elements_new[elements>cumorder[index]] - database[index,4] + database[index,3] - 1
      }
    }

    elements <- elements_new
    database <- database[which(activ_database),,drop=FALSE]

    if(TRUE){
      cumorder <- cumsum(c(1,diff(t(database[,3:4, drop=FALSE]))+1))
      elements_new <- elements
      for(index in 1:nrow(database)){
        remain <- intersect(elements, cumorder[index]:(cumorder[index+1]-1)) - cumorder[index] +1

        if(max(remain) < (database[index,4] - database[index,3]+1)){
          elements_new[elements>(cumorder[index+1]-1)] <- elements_new[elements>(cumorder[index+1]-1)] + max(remain) - (database[index,4] - database[index,3] +1 )
          database[index,4] <- max(remain) + database[index,3] - 1

        }
      }
      elements <- elements_new
    }
    elements <- elements_new

  } else{

    if(elements.copy){
      elements = which(!duplicated(get.id(population, database = database)))
    } else{
      elements <- 1:sum(database[,4]-database[,3]+1)
    }

  }
  if(depth.pedigree==Inf){
    pedigree.database <- get.database(population, gen=1:max(database[,1]))
  } else{
    new.pedigree.database <- new.pedigree.database_prior <- pedigree.database <- database
    remaining.depth <- depth.pedigree
    while(remaining.depth>0){
      parents <- get.pedigree(population, database = new.pedigree.database, raw=TRUE, include.error = include.error)
      m_parents <- rbind(parents[parents[,5]==1,4:6], parents[parents[,8]==1,7:9])
      f_parents <- rbind(parents[parents[,5]==2,4:6], parents[parents[,8]==2,7:9])
      if(nrow(m_parents)>0){
        m_gen <- unique(m_parents[,1])
        m_data <- cbind(m_gen, 1, 0,0)
        nincluded <- numeric(length(m_gen))
        for(index in 1:length(m_gen)){
          m_data[index,3] <- min(m_parents[m_parents[,1]==m_gen[index],3])
          m_data[index,4] <- max(m_parents[m_parents[,1]==m_gen[index],3])
          nincluded[index] <- length(unique(m_parents[m_parents[,1]==m_gen[index],3]))
        }

        for(index in length(m_gen):1){
          if(nincluded[index] < (m_data[index,4]-m_data[index,3]+1)/storage.save){
            m_data <- m_data[-index,]
            activ_p <- unique(m_parents[m_parents[,1]==m_gen[index],3])
            m_data <- rbind(m_data, cbind(m_gen[index], 1, activ_p, activ_p))
          }
        }

      } else{
        m_data <- NULL
      }
      if(nrow(f_parents)>0){
        f_gen <- unique(f_parents[,1])
        f_data <- cbind(f_gen, 2, 0,0)
        nincluded <- numeric(length(f_gen))
        for(index in 1:length(f_gen)){
          f_data[index,3] <- min(f_parents[f_parents[,1]==f_gen[index],3])
          f_data[index,4] <- max(f_parents[f_parents[,1]==f_gen[index],3])
          nincluded[index] <- length(unique(f_parents[f_parents[,1]==f_gen[index],3]))
        }

        for(index in length(f_gen):1){
          if(nincluded[index] < (f_data[index,4]-f_data[index,3]+1)/storage.save){
            f_data <- f_data[-index,]
            activ_p <- unique(f_parents[f_parents[,1]==f_gen[index],3])
            f_data <- rbind(f_data, cbind(f_gen[index], 2, activ_p, activ_p))
          }
        }

      } else{
        f_data <- NULL
      }

      new.pedigree.database <- get.database(population, database=rbind(m_data,f_data))
      new.pedigree.database <- unique(new.pedigree.database)

      #dups = new.pedigree.database[,1] * 1000000000 + new.pedigree.database[,2] * 100000000 + new.pedigree.database[,3] * 10000 + new.pedigree.database[,4]
      #dups2 = pedigree.database[,1]    * 1000000000 + pedigree.database[,2] * 100000000 + pedigree.database[,3] * 10000 + pedigree.database[,4]
      #already = dups %in% dups2
      #new.pedigree.database = new.pedigree.database[!already,,drop = FALSE]

      remaining.depth <- remaining.depth - 1
      pedigree.database <- rbind(new.pedigree.database, pedigree.database)

      if(length(new.pedigree.database) == 0 || (length(new.pedigree.database) == length(new.pedigree.database_prior) && prod(new.pedigree.database == new.pedigree.database_prior)==1)){
        break
      }
      new.pedigree.database_prior = new.pedigree.database
    }

    pedigree.database <- get.database(population, database = unique(pedigree.database))
  }

  ids_database <- get.id(population, database = database)
  ids_database_unique <- unique(ids_database)

  ids_pedigree = get.id(population, database = pedigree.database)
  if(min(diff(ids_pedigree))<1){
    ids_pedigree <- sort(unique(ids_pedigree))
  }
  ids_pedigree_first <- max(get.id(population, database = pedigree.database[pedigree.database[1,1]==pedigree.database[,1],,drop=FALSE]))

  n.animals <- length(ids_database_unique)
  n.total <- length(ids_pedigree)



  if(max(ids_pedigree) < 10^9 ){
    translate2 = numeric(max(ids_pedigree))

    translate2[ids_pedigree] = 1:length(ids_pedigree)

    position.pedigree = translate2[ids_database]
  } else{
    position.pedigree <- numeric(n.animals)
    for(index in 1:length(ids_database)){
      position.pedigree[index] <- which(ids_pedigree==ids_database[index])
    }
  }



  if(verbose) cat("Derive pedigree-matrix based for ", n.animals, " individuals based on ", n.total, " individuals.\n")
  kinship <- matrix(0L, ncol=n.total, nrow=n.total)

  group.size <- pedigree.database[,4]-pedigree.database[,3] +1
  if(length(start.kinship)==0){
    size.firstgen <- which(ids_pedigree_first==ids_pedigree)
    if(length(intersect(pedigree.database[1,1], population$info$founder.kinship))>0){

      activ_gen <- pedigree.database[1,1]
      if(ncol(population$info$kinship[[activ_gen]]) != sum(population$info$size[activ_gen,])){
        stop("Dimension of start kinship matrix does not match with generation size!")
      }
      founder_id <- get.id(population, gen=activ_gen)

      keeps = which(founder_id %in% ids_pedigree[1:size.firstgen])

      temp_kinship <- population$info$kinship[[activ_gen]][keeps, keeps] * int_mult2
      storage.mode(temp_kinship) <- "integer"
      kinship[1:size.firstgen, 1:size.firstgen] <- temp_kinship
    } else{
      kinship[1:size.firstgen, 1:size.firstgen] <- diag(as.integer(1/2 * int_mult),size.firstgen)
    }

  } else{
    kinship[1:nrow(start.kinship), 1:nrow(start.kinship)] <- start.kinship
    size.firstgen <- nrow(start.kinship)
  }
  first_new <- size.firstgen +1

  ## Potential export individual id in the pedigree - more efficient for high number of copies!
  #info.indi <- get.pedigree(population, database=pedigree.database)
  info.indi_id <- get.pedigree(population, database=pedigree.database, id=TRUE, include.error = include.error)

  #info.indi <-  info.indi[!duplicated(info.indi_id[,1]),,drop = FALSE]
  info.indi_id <-  info.indi_id[!duplicated(info.indi_id[,1]),,drop = FALSE]
  info.indi.pos <- matrix(0, ncol=3, nrow=nrow(info.indi_id))

  storage.mode(info.indi_id) = "numeric"

  if(max(ids_pedigree) < 10^9 ){
    translate = numeric(max(ids_pedigree))
    translate[ids_pedigree] = 1:length(ids_pedigree)
    info.indi.pos <- matrix(0, ncol=3, nrow=nrow(info.indi_id))
    info.indi.pos[info.indi_id!=0] = translate[info.indi_id]
  } else{
    for(index in 1:length(ids_pedigree)){
      info.indi.pos[info.indi_id==ids_pedigree[index]] <- index
    }
  }


  sorting <- sort(info.indi.pos[,1], index.return=TRUE)$ix

  #info.indi <- info.indi[sorting,,drop = FALSE]
  info.indi_id <- info.indi_id[sorting,,drop = FALSE]
  info.indi.pos <- info.indi.pos[sorting,,drop = FALSE]

  nr_indi <- info.indi.pos[,1]
  nr_father <- info.indi.pos[,2]
  nr_mother <- info.indi.pos[,3]


  if(first_new<= n.total){
    last <- 0
    for(second in first_new:n.total){
      nr.father <- nr_father[second]
      nr.mother <- nr_mother[second]
      first <- 1:second
      if(nr.father==0 && nr.mother==0){
        kinship[second,second] <- int_mult2
      } else{
        kinship[first,second] <- kinship[second,first] <- as.integer(0.5 * (if(nr.father==0){0} else{kinship[first, nr.father]} + if(nr.mother==0){0} else{kinship[first, nr.mother]}))
      }

      if(nr.mother==0 || nr.father==0 || (nr.father==nr.mother && nr.father==nr_indi[second])){
        kinship[second,second] <- int_mult2
        # Individual is founder!
      } else{
        kinship[second,second] <- int_mult2 + as.integer(0.5 * kinship[nr.father, nr.mother])
      }
    }
  }


  if(length(mult)>0){
    kinship.relevant <- kinship[position.pedigree,position.pedigree,drop = FALSE] / (int_mult /mult)
  } else{
    kinship.relevant <- kinship[position.pedigree,position.pedigree,drop = FALSE] / int_mult
  }

  kinship.relevant <- kinship.relevant[elements,elements,drop = FALSE]

  colnames(kinship.relevant) <- rownames(kinship.relevant) <- info.indi_id[position.pedigree[elements],1]

  t2 <- as.numeric(Sys.time())

  if(verbose){paste0(round(t2-t1, digits=2), " seconds for calculation of A matrix\n")}
  return(kinship.relevant)

}


#' Derive pedigree relationship matrix
#'
#' Function to derive pedigree matrix
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param depth.pedigree Depth of the pedigree in generations
#' @param start.kinship Relationship matrix of the individuals in the first considered generation
#' @param elements Vector of individuals from the database to include in pedigree matrix
#' @param storage.save Lower numbers will lead to less memory but slightly higher computing time (default: 1.5, min: 1)
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' data(ex_pop)
#' pedigree_matrix <- pedigree.matrix(population=ex_pop, gen=2)
#' @return Pedigree-based kinship matrix for in gen/database/cohort selected individuals
#' @export
#'

pedigree.matrix <- function(population, gen=NULL, database=NULL, cohorts=NULL, depth.pedigree=7,
                        start.kinship=NULL,
                        elements = NULL,
                        storage.save=1.5,
                        verbose=TRUE){

  A <- kinship.exp(population = population, gen=gen, database=database,
                   cohorts=cohorts, depth.pedigree=depth.pedigree,
                   start.kinship=start.kinship,
                   elements=elements,
                   storage.save=storage.save,
                   verbose=verbose,
                   mult=2)

  return(A)
}
