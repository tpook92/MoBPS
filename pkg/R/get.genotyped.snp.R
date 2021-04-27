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

#' Derive which markers are genotyped of selected individuals
#'
#' Function to devide which markers are genotyped for the selected individuals
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param export.alleles If TRUE export underlying alleles instead of just 012
#' @examples
#' data(ex_pop)
#' geno <- get.genotyped.snp(ex_pop, gen=2)
#' @return Binary Coded is/isnot genotyped level for in gen/database/cohorts selected individuals
#' @export

get.genotyped.snp <- function(population, database=NULL, gen=NULL, cohorts=NULL, export.alleles=FALSE){

  database <- get.database(population, gen, database, cohorts)


  start.chromo <- cumsum(c(1,population$info$snp)[-(length(population$info$snp)+1)])
  end.chromo <- population$info$cumsnp

  relevant.snps <- NULL
  chromosomen <- 1:length(population$info$snp)
  for(index in chromosomen){
    relevant.snps <- c(relevant.snps, start.chromo[index]:end.chromo[index])
  }
  nsnp <- length(relevant.snps)

  titel <- t(population$info$snp.base[,relevant.snps])

  n.animals <- sum(database[,4] - database[,3] +1)
  data <- matrix(0, ncol=n.animals, nrow=nsnp)
  before <- 0
  names <- numeric(n.animals)
  colnames(titel) <- c("Major_Allel", "Minor_Allel")

  for(row in 1:nrow(database)){
    animals <- database[row,]
    nanimals <- database[row,4] - database[row,3] +1

    if(nanimals>0){

        rindex <- 1
        for(index in animals[3]:animals[4]){
          arrays_used <- population$breeding[[animals[1]]][[animals[2]]][[index]][[22]]

          marker <- rep(FALSE, nsnp)
          for(ccc in arrays_used){
            marker[which(population$info$array.markers[[ccc]]==1)] <- TRUE
          }
          data[, before + rindex] <- marker
          rindex <- rindex + 1
        }
        before <- before + nanimals
      }
    }

  colnames(data) <- names
  rownames(data) <- population$info$snp.name[relevant.snps]
  if(export.alleles){
    return(list(titel,data))
  } else{
    return(data)
  }
}
