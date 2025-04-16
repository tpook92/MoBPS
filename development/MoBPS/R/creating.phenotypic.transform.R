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

#' Create a phenotypic transformation
#'
#' Function to perform create a transformation of phenotypes
#' @param population Population list
#' @param phenotypic.transform.function Phenotypic transformation to apply
#' @param trait Trait for which a transformation is to be applied
#' @param test.h2 Set to FALSE to not perform heritability check
#' @param h2 Vector of heritability input to test (before introducing noise from trafo; default: seq(0.05,0.5, by = 0.05))
#' @param export.h2 Set TRUE to export matrix of heritability before/after transformation
#' @param n.sample Sample size to use in test.h2 (default: 1000)
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' data(ex_pop)
#' trafo <- function(x){return(x^2)}
#' ex_pop <- creating.phenotypic.transform(ex_pop, phenotypic.transform.function=trafo)
#' @return Population-list with a new phenotypic transformation function
#' @export

creating.phenotypic.transform <- function(population, phenotypic.transform.function=NULL, trait=1,
                                          test.h2 = TRUE,
                                          gen = NULL,
                                          database = NULL,
                                          cohorts = NULL,
                                          h2 = seq(0.05,0.5, by = 0.05),
                                          export.h2 = FALSE,
                                          n.sample = 1000){

  if(length(phenotypic.transform.function)>0){
    population$info$phenotypic.transform[trait] <- TRUE
    population$info$phenotypic.transform.function[[trait]] <- phenotypic.transform.function
  }

  if(test.h2 || export.h2){

    if(length(h2)==0){
      h2 = 0
    }
    h2_obs = rep(0, length(h2))
    for(index in 1:length(h2)){
      n = 0

      pheno = NULL
      bv = NULL

      n_pheno = rep(0, population$info$bv.nr)
      n_pheno[trait] = 1

      if(length(gen)==0 && length(database)==0 && length(cohorts) == 0){
        gen = get.ngen(population)
      }
      database = get.database(population, gen = gen, database = database, cohorts = cohorts)

      database[  database[,4] > 1000,4] = n.sample
      while(n < n.sample){

        tmp3 = h2[index]
        if(tmp3 == 0){
          tmp3 = NULL
        }
        pop1 = breeding.diploid(population, n.observation = n_pheno,
                                heritability = tmp3,
                                phenotyping.database = database,
                                verbose = FALSE)

        pheno = c(pheno, get.pheno(pop1, database = database)[trait,])
        bv = c(bv, get.bv(pop1, database = database)[trait,])
        n = length(bv)
      }
      h2_obs[index] = stats::cor(pheno, bv)^2
    }

    cat("Heritabilities provided downstream will be for underlying quantitative trait (not discrete phenotypes)\n")
    cat("Expected heritabilities:\n")
    tmp = cbind(h2, h2_obs)
    colnames(tmp) = c("Underlying Quantitative Trait", "Observed Phenotype")


    if(export.h2){
      return(tmp)
    } else{
      print(tmp)
    }

  }
  return(population)

}
