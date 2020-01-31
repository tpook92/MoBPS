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

#' Export location of individuals from the population list
#'
#' Export location of individuals from the population list
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @examples
#' data(ex_pop)
#' get.individual.loc(ex_pop, gen=2)
#' @return Storage Position for in gen/database/cohorts selected individuals (Generation/Sex/IndividualNr)
#' @export

get.individual.loc <- function(population, database=NULL, gen=NULL, cohorts=NULL){

  database <- get.database(population, gen, database, cohorts)

  n.animals <- sum(database[,4] - database[,3] +1)
  data <- matrix(0, nrow=n.animals, ncol=3)
  before <- 0
  rown <- numeric(n.animals)
  for(row in 1:nrow(database)){
    animals <- database[row,]
    nanimals <- database[row,4] - database[row,3] +1
    data[(before+1):(before+nanimals),] <- cbind(database[row,1], database[row,2], database[row,3]:database[row,4])
    rown[(before+1):(before+nanimals)] <- paste(if(database[row,2]==1) "M" else "F", database[row,3]:database[row,4], "_", database[row,1], sep="")
    before <- before + nanimals

  }
  colnames(data) <- c("generation", "sex", "individual nr.")
  row.names(data) <- rown

  return(data)

}
