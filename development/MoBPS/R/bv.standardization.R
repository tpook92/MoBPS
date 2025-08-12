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

#' BV standardization
#'
#' Function to get mean and genetic variance of a trait to a fixed value
#' @param population Population list
#' @param mean.target Target mean
#' @param var.target Target variance
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param adapt.bve Modify previous breeding value estimations by scaling (default: FALSE)
#' @param adapt.pheno Modify previous phenotypes by scaling (default: FALSE)
#' @param adapt.sigma.e Set to TRUE to scale sigma.e values used based on scaling
#' @param verbose Set to TRUE to display prints
#' @param set.zero Set to TRUE to have no effect on the 0 genotype (or 00 for QTLs with 2 underlying SNPs)
#' @param traits Use this parameter to only perform scaling of these traits (alternatively set values in mean/var.target to NA, default: all traits)
#' @examples
#' population <- creating.diploid(nsnp=1000, nindi=100, n.additive=100)
#' population <- bv.standardization(population, mean.target=200, var.target=5)
#' @return Population-list with scaled QTL-effects
#' @export



bv.standardization <- function(population, mean.target=NA, var.target=NA, gen=NULL, database=NULL, cohorts=NULL,
                               adapt.bve=TRUE, adapt.pheno=NULL, verbose=FALSE, set.zero = FALSE,
                               adapt.sigma.e = FALSE,
                               traits = NULL){



  n_traits <- population$info$bv.nr

  if(length(traits)>0){
    mean.target_temp = rep(NA, n_traits)
    var.target_temp = rep(NA, n_traits)
    mean.target_temp[traits] = mean.target
    var.target_temp[traits] = var.target
    mean.target = mean.target_temp
    var.target = var.target_temp
  }

  if(length(adapt.pheno)==0){
    if(sum(population$info$phenotypic.transform)>0){
      adapt.pheno = FALSE
      if(verbose){
        cat("Phenotype transformation deactivated as phenotypic transformation is used. Set adapt.pheno = TRUE to scale.\n")
      }
    } else{
      adapt.pheno = TRUE
    }
  }

  modi1 <- rep(1, n_traits)
  modi2 <- population$info$base.bv

  if(length(mean.target)<n_traits){
    mean.target <- rep(mean.target, length.out = n_traits)
    mean.target <- c(mean.target, rep(NA, length(mean.target)-n_traits))
  }
  if(length(var.target)<n_traits){
    var.target <- rep(var.target, length.out = n_traits)
    var.target <- c(var.target, rep(NA, length(var.target)-n_traits))
  }
  if(length(gen)==0 && length(database)==0 && length(cohorts)==0){
    gen <- nrow(population$info$size)
  }
  database <- get.database(population, gen, database, cohorts)
  ## Variance Standardization
  for(index in (1:n_traits)[!population$info$is.combi]){
    if(!is.na(var.target[index])){
      new_var <- var.target[index]

      if(population$info$bv.calculated==FALSE){
        population <- breeding.diploid(population, verbose=verbose)
      }

      var_test <- stats::var(get.bv(population, database= database)[index,])

      if(var_test==0){
        stop(paste0("No variance in trait ", index, ". No scaling to achieve target variance possible."))
      }

      test1 <- TRUE
      if(length(population$info$real.bv.add[[index]])>0){
        population$info$real.bv.add[[index]][,3:5] <- population$info$real.bv.add[[index]][,3:5] * sqrt(  new_var / var_test)
        if(set.zero){
          population$info$real.bv.add[[index]][,3:5] = population$info$real.bv.add[[index]][,3:5] - population$info$real.bv.add[[index]][,3]
        }
        test1 <- FALSE
      }
      if(length(population$info$real.bv.mult[[index]])>0){
        population$info$real.bv.mult[[index]][,5:13] <- population$info$real.bv.mult[[index]][,5:13] * sqrt(  new_var / var_test)
        if(set.zero){
          population$info$real.bv.mult[[index]][,5:13] = population$info$real.bv.mult[[index]][,5:13] - population$info$real.bv.mult[[index]][,5]
        }
        test1 <- FALSE
      }
      modi1[index] <- sqrt(new_var / var_test)

      if(test1 && verbose) cat("You entered a trait without quantitative loci. Is this intentional?\n")
    }


  }

  for(gen in 1:length(population$breeding)){
    for(sex in 1:2){
      if(length(population$breeding[[gen]][[sex+6]])>0){
        population$breeding[[gen]][[sex+6]] <- (((population$breeding[[gen]][[sex+6]] - modi2) * modi1) + modi2)
      }
    }
  }


  population$info$pool_effects_calc = FALSE
  population$info$pool_list = NULL
  population$info$bypool_list = NULL

  ## Mean Standardization
  for(index in 1:n_traits){

    if(!is.na(mean.target[index])){
      if(population$info$bv.calculated==FALSE){
        population <- breeding.diploid(population, verbose=FALSE)
      }

      mean_test <- mean(get.bv(population, database = database)[index,])

      population$info$base.bv[index] <- mean.target[index] + population$info$base.bv[index] - mean_test
    }


  }

  for(gen in 1:length(population$breeding)){
    for(sex in 1:2){
      if(length(population$breeding[[gen]][[sex+6]])>0){
        population$breeding[[gen]][[sex+6]] <- population$breeding[[gen]][[sex+6]] - modi2 + population$info$base.bv
      }
    }
  }


  if(adapt.bve){
    for(gen in 1:length(population$breeding)){
      for(sex in 1:2){
        activ <- (population$breeding[[gen]][[sex+2]]!= 0)
        if(sum(activ)>0){
          population$breeding[[gen]][[sex+2]][activ] <- (((population$breeding[[gen]][[sex+2]] - modi2) * modi1) + population$info$base.bv)[activ]
        }
      }
    }
  }
  if(adapt.pheno){
    for(gen in 1:length(population$breeding)){
      for(sex in 1:2){
        activ <- !is.na(population$breeding[[gen]][[sex+8]]) & (population$breeding[[gen]][[sex+8]]!= 0) & (!population$info$phenotypic.transform)

        if(sum(activ)>0){
          population$breeding[[gen]][[sex+8]][activ] <- (((population$breeding[[gen]][[sex+8]] - modi2) * modi1) + population$info$base.bv)[activ]
        }
      }
    }
  }

  if(adapt.sigma.e){
    population$info$last.sigma.e.value = population$info$last.sigma.e.value * modi1
  }

  if(length(population$info$e0)>0){
    population$info$e0 = NULL
  }
  if(length(population$info$e1)>0){
    population$info$e1 = NULL
  }
  if(length(population$info$e2)>0){
    population$info$e2 = NULL
  }


  return(population)
}
