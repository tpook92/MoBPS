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

#' Generation of genomic traits
#'
#' Generation of the trait in a starting population
#' @param population Population list
#' @param merge Vector containing the traits to merge (e.g. c(2,3))
#' @param trait.name Name of the trait generated
#' @param bv.standard Set TRUE to standardize trait mean and variance via bv.standardization()
#' @param mean.target Target mean
#' @param var.target Target variance
#' @param verbose Set to FALSE to not display any prints
#' @param set.zero Set to TRUE to have no effect on the 0 genotype (or 00 for QTLs with 2 underlying SNPs)
#' @param new.phenotype.correlation (OLD! - use new.residual.correlation) Correlation of the simulated enviromental variance
#' @param new.residual.correlation Correlation of the simulated enviromental variance
#' @examples
#' population <- creating.diploid(nsnp=1000, nindi=100)
#' population <- creating.trait(population, n.additive=c(100,100))
#' population <- merging.trait(population, merge = c(1,2))
#' @return Population-list with one or more additional new traits
#' @export


merging.trait <- function(population,
                        merge = NULL,
                           trait.name=NULL,
                           bv.standard=FALSE,
                           mean.target=NULL,
                           var.target=NULL,
                           verbose=TRUE,
                           set.zero = FALSE,
                           new.phenotype.correlation = NULL,
                        new.residual.correlation = NULL){

  if(length(merge)==0){
    stop("Please provide input for which traits to merge in merge parameter")
  }

  if(length(mean.target)>0){
    bv.standard <- TRUE
  } else{
    mean.target <- 100
  }
  if(length(var.target)>0){
    bv.standard <- TRUE
  } else{
    var.target <- 10
  }

  if(sum(set.zero)>0){
    bv.standard <- TRUE
  }

  keep = min(merge)
  delete = merge[merge!=keep]

  ntrait_before = population$info$bv.nr
  ntrait = ntrait_before - length(delete)


  preserve.bve <- length(population)==0

  if(length(new.phenotype.correlation)>0){
    new.residual.correlation <- new.phenotype.correlation
  }


  real.bv.add_temp = NULL
  real.bv.mult_temp = NULL
  real.bv.dice_temp = list(list(), list())

  for(bven in merge){
    if(length(population$info$real.bv.add[[bven]])>0){
      real.bv.add_temp = rbind(real.bv.add_temp, population$info$real.bv.add[[bven]])
    }

    if(length(population$info$real.bv.mult[[bven]])>0){
      real.bv.mult_temp = rbind(real.bv.mult_temp, population$info$real.bv.mult[[bven]])
    }
    if(length(population$info$real.bv.dice[[bven]][[1]])>0){
      real.bv.dice_temp[[1]] == c(real.bv.dice_temp[[1]] , population$info$real.bv.dice[[bven]][[1]])
      real.bv.dice_temp[[2]] == c(real.bv.dice_temp[[2]] , population$info$real.bv.dice[[bven]][[2]])
    }

  }


  if(length(real.bv.add_temp)>0){
    population$info$real.bv.add[[keep]] = real.bv.add_temp
  }
  if(length(real.bv.mult_temp)>0){
    population$info$real.bv.mult[[keep]] = real.bv.mult_temp
  }
  if(length(real.bv.dice_temp[[1]])>0){
    population$info$real.bv.dice[[keep]] = real.bv.dice_temp
  }

  for(bven in delete){
    population$info$real.bv.add[[bven]] = NULL
    population$info$real.bv.mult[[bven]] = NULL
    population$info$real.bv.dice[[bven]] = NULL
  }

  population$info$bv.nr =  population$info$bv.nr - length(delete)
  population$info$bv.calc = population$info$bv.calc - sum(delete <=population$info$bv.calc )

  # this is not all necessary but does not hurt
  if(TRUE){
    population$info$bv.calculated.partly <- NULL
  }

  bv.random = population$info$bv.random
  bv.random.variance = population$info$bv.random.variance
  bv.random.variance[keep] = sum(bv.random.variance[merge])

  if(length(unique(bv.random))>1){
    stop("Merging of QTL traits with non-QTL based traits is not allowed")
  }

  bv.random = bv.random[-delete]
  bv.random.variance = bv.random.variance[-delete]

  population$info$bv.random <- bv.random
  population$info$bv.random.variance <- bv.random.variance


  if(population$info$phenotypic.transform[keep]){
    warning(paste0("Use phenotypic transformation function of trait ", keep))
  }
  population$info$phenotypic.transform <- population$info$phenotypic.transform[-delete]
  for(bven in delete){
    population$info$phenotypic.transform.function[[bven]] <- NULL
  }

  if(length(unique(population$info$is.maternal))>1){
    stop("Merging of materal trait with regular trait is not allowed!")
  }
  population$info$is.maternal = population$info$is.maternal[-delete]

  if(length(unique(population$info$is.paternal))>1){
    stop("Merging of pateral trait with regular trait is not allowed!")
  }
  population$info$is.paternal = population$info$is.paternal[-delete]

  if(sum(population$info$is.combi)!=0){
    stop("Merging of combi traits is not allowed!")
  }
  population$info$is.combi = population$info$is.combi[-delete]


  if(length(unique(population$info$bve.mult.factor))>1){
    warning(paste0("Use multiplication factor of trait ", keep))
  }
  population$info$bve.mult.factor = population$info$bve.mult.factor[-delete]

  if(length(unique(population$info$bve.poly.factor))>1){
    warning(paste0("Use polynomical scaling of trait ", keep))
  }
  population$info$bve.poly.factor = population$info$bve.poly.factor[-delete]

  population$info$base.bv[keep] = sum( population$info$base.bv[merge])
  population$info$base.bv = population$info$base.bv[-delete]


  if(length(new.residual.correlation)==0){
    P_new = diag(ntrait_before)
    P_new[merge,keep] <- 1
    P_new[,delete] = 0

    new.residual.correlation = (t(P_new) %*% population$info$pheno.correlation) %*% (t(population$info$pheno.correlation) %*% (P_new))

    new.residual.correlation = diag(1/sqrt(diag(new.residual.correlation))) %*% new.residual.correlation %*% diag(1/sqrt(diag(new.residual.correlation)))

    new.residual.correlation = new.residual.correlation[-delete, -delete, drop = FALSE]
  }



  new.breeding.correlation = (t(P_new) %*% population$info$bv.correlation) %*% ((P_new))

  new.breeding.correlation = diag(1/sqrt(diag(new.breeding.correlation))) %*% new.breeding.correlation %*% diag(1/sqrt(diag(new.breeding.correlation)))

  new.breeding.correlation = new.breeding.correlation[-delete, -delete, drop = FALSE]

  population$info$bv.correlation <- new.breeding.correlation

  population$info$trait.name = population$info$trait.name[-delete]
  if(length(trait.name)==1){
    population$info$trait.name[keep] == trait.name
  }

  for(generation in 1:nrow(population$info$size)){

    for(tt in c(3,4,7:10,19:22,27:30)){
      population$breeding[[generation]][[tt]][keep,] <- 0
      population$breeding[[generation]][[tt]] <- population$breeding[[generation]][[tt]][-delete,,drop=FALSE]
    }
  }


  population$info$fixed.effects[keep,] = colSums(population$info$fixed.effects[merge,,drop = FALSE])
  population$info$fixed.effects = population$info$fixed.effects[-delete,,drop = FALSE]



  # Add traits with no generated phenotypes / litter or pen effects.
  temp1 <- rep(0, population$info$bv.nr)
  for(gen in 1:length(population$breeding)){
    for(sex in 1:2){
      if(length(population$breeding[[gen]][[sex]])>0){
        for(index in 1:length(population$breeding[[gen]][[sex]])){
          population$breeding[[gen]][[sex]][[index]][[15]] <- temp1
          population$breeding[[gen]][[sex]][[index]][[28]] <- c(population$breeding[[gen]][[sex]][[index]][[28]], rep(0, ncol(population$info$fixed.effects) - length(population$breeding[[gen]][[sex]][[index]][[28]])))
          population$breeding[[gen]][[sex]][[index]][[29]] <- temp1
          population$breeding[[gen]][[sex]][[index]][[30]] <- temp1
        }
      }
    }
  }

  population$info$real.bv.length[1] <- length(population$info$real.bv.add) - 1
  population$info$real.bv.length[2] <- length(population$info$real.bv.mult) - 1
  population$info$real.bv.length[3] <- length(population$info$real.bv.dice) - 1

  if(bv.standard){
    population <- bv.standardization(population, mean.target = mean.target, var.target = var.target, set.zero = set.zero)
  }


  # Calculation of initial genomic values
  population$info$bv.calculated = FALSE
  population <- breeding.diploid(population, verbose=FALSE)

  return(population)
}
