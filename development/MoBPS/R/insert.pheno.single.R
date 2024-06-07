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

#' Manually enter estimated breeding values
#'
#' Function to manually enter estimated breeding values
#' @param population Population list
#' @param phenos Matrix of breeding values to enter (one row per individual: 1st column: individual, 2nd column: trait nr, 3rd-nth column records)
#' @param count.only.increase Set to FALSE to reduce the number of observation for a phenotype to "count" (default: TRUE)
#' @examples
#' data(ex_pop)
#' bv <- get.bv(ex_pop, gen=2)
#' new.bve <- cbind( colnames(bv), bv[,1]) ## Unrealistic but you do not get better than this!
#' ex_pop <- insert.bve(ex_pop, bves=new.bve)
#' @return Population-List with newly entered estimated breeding values
#' @export

insert.pheno.single <- function(population, phenos, count.only.increase=TRUE){

  bves = phenos

  add <- 8

  sex <- as.numeric(substr(bves[,1], start=1, stop=1)=="F")+1
  split <- unlist(strsplit(bves[,1], split=c("_")))
  gen <- as.numeric(split[1:length(sex)*2])
  nr <- as.numeric(substr(split[1:length(sex)*2-1], start=2, stop=100))

  databases <- cbind(gen, sex, nr)

  groups <- unique(databases[,1:2])

  for(index in 1:nrow(groups)){
    activ <- databases[,1]==groups[index,1] & databases[,2]==groups[index,2]
    database_activ <- databases[activ,,drop = FALSE]
    bves_activ <- t(bves[activ,-1, drop = FALSE])
    storage.mode(bves_activ) <- "numeric"



    for(index2 in 1:nrow(database_activ)){
      sex <- groups[index,2]
      gen <- groups[index,1]
      nr <- database_activ[index2,3]
      bven = bves_activ[1,index2]
      tmp1 = bves_activ[-1,index2]

      if(!count.only.increase || sum(!is.na(tmp1)) >= population$breeding[[gen]][[sex]][[nr]][[15]][bven]){
        population$breeding[[gen]][[sex]][[nr]][[15]][bven] = sum(!is.na(tmp1))
        population$breeding[[gen]][[sex+add]][bven,nr] <- mean(tmp1, na.rm = TRUE)

        population$breeding[[gen]][[sex]][[nr]][[27]][[bven]] = tmp1[!is.na(tmp1)]

        max_obs =  max(population$breeding[[gen]][[sex]][[nr]][[15]])
        if(length(population$breeding[[gen]][[sex]][[nr]][[24]])>0 && ncol(population$breeding[[gen]][[sex]][[nr]][[24]])>max_obs){
          if(max_obs>0){
            population$breeding[[gen]][[sex]][[nr]][[24]] <- population$breeding[[gen]][[sex]][[nr]][[24]][,1:max_obs]
          } else{
            population$breeding[[gen]][[sex]][[nr]][24] <- list(NULL) ## Only single bracket to not reduce length of list
          }

        }
      }


    }


  }

  return(population)
}
