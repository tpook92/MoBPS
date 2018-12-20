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

#' Bestimmung der Urspruenge ausgewaehlter Tiere
#'
#' Funktion zur Bestimmung der Ursprungstiere fuer das genetische Material ausgewaehlter Tiere.
#' @param population Datenvektor
#' @param database Menge ausgewahlter Tiere - Matrix( pro Spalte: Generation / Geschlecht)
#' @param gen Schnelleingabe von database (Vektor mit allen relevanten Generationen)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @export

get.recombi <- function(population, database=NULL, gen=NULL, cohorts=NULL){

  if(length(gen)>0){
    database <- cbind(rep(gen,each=2), rep(1:2, length(gen)))
  }
  if(length(database)>0 && ncol(database)==2){
    start <- end <- numeric(nrow(database))
    for(index in 1:nrow(database)){
      start[index] <- 1
      end[index] <- population$info$size[database[index,1], database[index,2]]
    }
    database <- cbind(database, start, end)
  }
  if(length(cohorts)>0){
    database2 <- matrix(0, nrow=length(cohorts), ncol=4)
    for(index in 1:length(cohorts)){
      row <- which(population$info$cohorts==cohorts[index])
      gen <- as.numeric(population$info$cohorts[row,2])
      sex <- 1 + (as.numeric(population$info$cohorts[row,4])>0)
      first <- as.numeric(population$info$cohorts[row,5 + sex])
      last <- first + as.numeric(population$info$cohorts[row,2 + sex]) - 1
      database2[index,] <- c(gen,sex,first,last)
    }
    database <- rbind(database, database2)
  }

  n.animals <- sum(database[,4] - database[,3] +1)

  data <- list()
  rindex <- 1

  for(row in 1:nrow(database)){
    animals <- database[row,]
    for(index in database[row,3]:database[row,4]){
      data[[rindex]] <-  list()

      colnamed <- paste(if(animals[2]==1) "M" else "W", index,"_", animals[1],sep="")
      data[[rindex]][[1]] <- population$breeding[[animals[1]]][[animals[2]]][[index]][[1]]
      data[[rindex]][[2]] <- population$breeding[[animals[1]]][[animals[2]]][[index]][[2]]
      if(population$info$miraculix){
        if (requireNamespace("miraculix", quietly = TRUE)) {
          decode <- miraculix::decodeOrigins
        } else{
          decode <- decodeOriginsR
        }
      } else{
        decode <- decodeOriginsR
      }
      if(population$info$miraculix){
        male_ori <- matrix(0, nrow=length(population$breeding[[animals[1]]][[animals[2]]][[index]][[5]]), ncol=4)
        female_ori <- matrix(0, nrow=length(population$breeding[[animals[1]]][[animals[2]]][[index]][[6]]), ncol=4)
        for(index2 in 1:nrow(male_ori)){
          male_ori[index2,] <- decode(population$breeding[[animals[1]]][[animals[2]]][[index]][[5]],index2)
        }
        for(index2 in 1:nrow(female_ori)){
          female_ori[index2,] <- decode(population$breeding[[animals[1]]][[animals[2]]][[index]][[6]],index2)
        }
        data[[rindex]][[3]] <- male_ori
        data[[rindex]][[4]] <- female_ori

      } else{
        data[[rindex]][[3]] <- matrix(decode(population$breeding[[animals[1]]][[animals[2]]][[index]][[5]]), byrow=FALSE, ncol=4)
        data[[rindex]][[4]] <- matrix(decode(population$breeding[[animals[1]]][[animals[2]]][[index]][[6]]), byrow=FALSE, ncol=4)

      }

      data[[rindex]][[5]] <- colnamed


      rindex <- rindex + 1
    }
  }


  return(data)
}
