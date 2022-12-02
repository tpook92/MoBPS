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

#' write.pedigree.mixblup
#'
#' write.pedigree.mixblup
#' @param population Population list
#' @param path AA
#' @param gen AA
#' @param database AA
#' @param cohorts AA
#' @param id AA
#' @param depth.pedigree AA
#' @param storage.save AA
#' @param verbose AA
#' @return write.pedigree.mixblup
#' @examples
#' write.pedigree.mixblup
#' @export


write.pedigree.mixblup <- function(population, path, gen=NULL, database=NULL, cohorts=NULL , id = NULL, depth.pedigree=7,
                           storage.save = 1.5, verbose=TRUE){

  if(verbose) cat(paste0("Start writting pedigree file at ", path,"\n"))
  database = get.database(population, gen = gen, database = database, cohorts = cohorts, id = id)

  if(verbose) cat(paste0("Start collecting pedigree\n"))
  # Generate data needed for the ped-file (see kinship.exp())
  if(depth.pedigree==Inf){
    pedigree.database <- get.database(population, gen=1:max(database[,1]))
  } else{
    new.pedigree.database <- pedigree.database <- database
    remaining.depth <- depth.pedigree
    while(remaining.depth>0){
      parents <- get.pedigree(population, database = new.pedigree.database, raw=TRUE)
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
      remaining.depth <- remaining.depth - 1
      pedigree.database <- rbind(new.pedigree.database, pedigree.database)
    }

    pedigree.database <- get.database(population, database = pedigree.database)
  }

  pedigree_table <- get.pedigree(population, database = pedigree.database, id=TRUE)

  add <- pedigree_table[which(!duplicated(as.character(pedigree_table))[-(1:nrow(pedigree_table))]) + nrow(pedigree_table)]
  add = add[add!=0]

  pedigree_table = pedigree_table[!duplicated(pedigree_table[,1]),]

  if(verbose) cat(paste0("Pedigree contains ", nrow(pedigree_table) + length(add), " animals with ", length(add), " animals without known parents\n"))
  if(verbose) cat(paste0("Start writting: ", path,"\n"))

  if (requireNamespace("fwrite", quietly = TRUE)) {
    if(length(add)>0){
      data.table::fwrite(file=path, cbind(add,0,0), col.names = FALSE)
      data.table::fwrite(file=path, pedigree_table, col.names = FALSE, append = TRUE)
    } else{
      data.table::fwrite(file=path, pedigree_table, col.names = FALSE)
    }
  } else{
    if(length(add)>0){
      utils::write.table(file=path, cbind(add,0,0), col.names = FALSE, row.names = FALSE, quote = FALSE)
      utils::write.table(file=path, pedigree_table, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
    } else{
      utils::write.table(file=path, pedigree_table, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
  }




}
