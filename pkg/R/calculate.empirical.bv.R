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

#' Export estimated breeding values
#'
#' Function to export estimated breeding values
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param id Individual IDs to search/collect in the database
#' @param breeding.size Number of offspring to generate for calculations (default: c(100,0))
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names (default: TRUE)
#' @param reference.database Groups of individuals to consider as mating partners
#' @param reference.gen Quick-insert for database (vector of all generations to export)
#' @param reference.cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param reference.id Individual IDs to search/collect in the database
#' @param reference.freq Per-SNP allelle frequencies of a reference population to mate to (default: NULL)
#' @param reference.dataset SNP dataset, use "random", "macs", "allhetero", "all0" to mate to (default: NULL)
#' @param enter.population Set to TRUE to enter empirical underlying true breeding values in the population list instead of exporting matrix
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' data(ex_pop)
#' calculate.empirical.bv(ex_pop, gen=1, reference.gen = 2)
#' @return Estimated breeding value of in gen/database/cohorts selected individuals
#' @export


calculate.empirical.bv = function(population, database=NULL, gen=NULL, cohorts=NULL, use.id=TRUE, id = NULL,
                 breeding.size = c(100,0), verbose = TRUE,
                 reference.gen = NULL, reference.database = NULL, reference.cohorts = NULL, reference.id = NULL,
                 reference.freq = NULL, reference.dataset = NULL,
                 enter.population = FALSE
                 ){

  database <- get.database(population, gen = gen, database = database,
                           cohorts = cohorts, id = id, per.individual = TRUE)

  if(length(reference.database) ==0 && length(reference.id)==0 && length(reference.gen)==0 && length(reference.cohorts)==0){

    population = creating.diploid(population, freq = reference.freq, dataset = reference.database,
                                  sex.quota = 1, verbose = FALSE,
                                  name.cohort = "MoPBS_Mating_pairs")

    reference.cohorts = "MoPBS_Mating_pairs"
  }

  reference.database <- get.database(population, gen = reference.gen, database = reference.database,
                           cohorts = reference.cohorts, id = reference.id)


  n.animals <- sum(database[,4] - database[,3] +1)
  data <- matrix(0, ncol=n.animals, nrow=population$info$bv.nr)

  if(verbose){
    options(scipen = 999)
    cat("Empirical calculation of true breeding values will require simulation of ", nrow(database) * sum(breeding.size), " matings.\n")
  }
  names <- numeric(n.animals)

  if(verbose){
    pb <- utils::txtProgressBar(min = 0, max = nrow(database), style = 3)
    pb_temp <- round(nrow(database)/25)
  }

  for(row in 1:nrow(database)){

    pop1 = breeding.diploid(population, breeding.size = c(100,0),
                            selection.m.database = database[row,],
                            selection.f.database = reference.database,
                            name.cohort = "MoBPS_TmpOffspring",
                            verbose = FALSE)

    data[,row] <- rowMeans(get.bv(pop1, cohorts = "MoBPS_TmpOffspring"))
    names[row] <- paste(if(database[row,2]==1) "M" else "F", database[row,3],"_", database[row,1], sep="")

    if(verbose & (nrow(database) < 25 || row%%pb_temp==0)){
      utils::setTxtProgressBar(pb, row)
    }

  }

  if(verbose ){
    close(pb)

  }

  rownames(data) <- population$info$trait.name

  if(use.id){
    colnames(data) <- get.id(population, database = database)
  } else{
    colnames(data) <- names
  }


  if(!enter.population){
    return(data)
  } else{

    for(row in 1:nrow(database)){


      db_tmp = mat3(population$breeding[[database[row,1]]][[database[row,2]]][[database[row,3]]][[21]])

      for(row2 in 1:nrow(db_tmp)){
        activ = db_tmp[row2,]

        if(length(population$breeding[[activ[1]]][[activ[2] + 48]]) == 0 ){

          population$breeding[[activ[1]]][[activ[2]+48]] = matrix(NA, nrow= population$info$bv.nr, ncol=population$info$size[activ[1], activ[2]])
        }

        population$breeding[[activ[1]]][[activ[2]+48]][,activ[3]] = data[,row]
      }



    }

    return(population)
  }

}
