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

#' Generation of genomic traits
#'
#' Generation of the trait in a starting population
#' @param population Population list
#' @param bv.total Number of traits (If more than traits via real.bv.X use traits with no directly underlying QTL)
#' @param polygenic.variance Genetic variance of traits with no underlying QTL
#' @param randomSeed Set random seed of the process
#' @param n.additive Number of additive QTL with effect size drawn from a gaussian distribution
#' @param n.dominant Number of dominant QTL with effect size drawn from a gaussian distribution
#' @param n.equal.additive Number of additive QTL with equal effect size (effect.size)
#' @param n.equal.dominant Number of n.equal.dominant QTL with equal effect size
#' @param n.qualitative Number of qualitative epistatic QTL
#' @param n.quantitative Number of quantitative epistatic QTL
#' @param dominant.only.positive Set to TRUE to always asign the heterozygous variant with the higher of the two homozygous effects (e.g. hybrid breeding); default: FALSE
#' @param var.additive.l Variance of additive QTL
#' @param var.dominant.l Variance of dominante QTL
#' @param var.qualitative.l Variance of qualitative epistatic QTL
#' @param var.quantitative.l Variance of quantitative epistatic QTL
#' @param effect.size.equal.add Effect size of the QTLs in n.equal.additive
#' @param effect.size.equal.dom Effect size of the QTLs in n.equal.dominant
#' @param exclude.snps Marker were no QTL are simulated on
#' @param replace.traits If TRUE delete the simulated traits added before
#' @param shuffle.traits Combine different traits into a joined trait
#' @param shuffle.cor Target Correlation between shuffeled traits
#' @param real.bv.add Single Marker effects
#' @param real.bv.mult Two Marker effects
#' @param real.bv.dice Multi-marker effects
#' @param bve.mult.factor Multiplicate trait value times this
#' @param bve.poly.factor Potency trait value over this
#' @param base.bv Average genetic value of a trait
#' @param new.phenotype.correlation (OLD! - use new.residual.correlation) Correlation of the simulated enviromental variance
#' @param new.residual.correlation Correlation of the simulated enviromental variance
#' @param new.breeding.correlation Correlation of the simulated genetic variance (child share! heritage is not influenced!
#' @param trait.name Name of the trait generated
#' @param remove.invalid.qtl Set to FALSE to deactive the automatic removal of QTLs on markers that do not exist
#' @param bv.standard Set TRUE to standardize trait mean and variance via bv.standardization()
#' @param mean.target Target mean
#' @param var.target Target variance
#' @param verbose Set to FALSE to not display any prints
#' @param is.maternal Vector coding if a trait is caused by a maternal effect (Default: all FALSE)
#' @param is.paternal Vector coding if a trait is caused by a paternal effect (Default: all FALSE)
#' @examples
#' population <- creating.diploid(nsnp=1000, nindi=100)
#' population <- creating.trait(population, n.additive=100)
#' @return Population-list with one or more additional new traits
#' @export


creating.trait <- function(population, real.bv.add=NULL, real.bv.mult=NULL, real.bv.dice=NULL,
                           bv.total=0, polygenic.variance=100,
                           bve.mult.factor=NULL, bve.poly.factor=NULL, base.bv=NULL,
                           new.phenotype.correlation=NULL,
                           new.residual.correlation=NULL,
                           new.breeding.correlation=NULL,
                           n.additive=0,
                           n.equal.additive=0,
                           n.dominant=0,
                           n.equal.dominant=0,
                           n.qualitative=0,
                           n.quantitative=0,
                           dominant.only.positive = FALSE,
                           var.additive.l=NULL,
                           var.dominant.l=NULL,
                           var.qualitative.l=NULL,
                           var.quantitative.l=NULL,
                           effect.size.equal.add = 1,
                           effect.size.equal.dom = 1,
                           exclude.snps=NULL,
                           randomSeed=NULL,
                           shuffle.traits=NULL,
                           shuffle.cor= NULL,
                           replace.traits=FALSE,
                           trait.name=NULL,
                           remove.invalid.qtl=TRUE,
                           bv.standard=FALSE,
                           mean.target=NULL,
                           var.target=NULL,
                           verbose=TRUE,
                           is.maternal=NULL,
                           is.paternal=NULL){


  if(length(randomSeed)>0){
    set.seed(randomSeed)
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


  preserve.bve <- length(population)==0

  if(length(new.phenotype.correlation)>0){
    new.residual.correlation <- new.phenotype.correlation
  }

  if(!is.list(var.additive.l) ){
    var.additive.l <- list(var.additive.l)
  }
  if(!is.list(var.dominant.l)){
    var.dominant.l <- list(var.dominant.l)
  }
  if(!is.list(var.qualitative.l)){
    var.qualitative.l <- list(var.qualitative.l)
  }
  if(!is.list(var.quantitative.l)){
    var.quantitative.l <- list(var.quantitative.l)
  }

  trait_sum <- n.additive + n.dominant + n.qualitative + n.quantitative
  test <- list(NULL)
  if(length(var.additive.l) < length(trait_sum)){
    var.additive.l <- c(var.additive.l, rep(test,length.out=length(trait_sum)-length(var.additive.l)))
  }
  if(length(var.dominant.l) < length(trait_sum)){
    var.dominant.l <- c(var.dominant.l, rep(test,length.out=length(trait_sum)-length(var.dominant.l)))
  }
  if(length(var.qualitative.l) < length(trait_sum)){
    var.qualitative.l <- c(var.qualitative.l, rep(test,length.out=length(trait_sum)-length(var.qualitative.l)))
  }
  if(length(var.quantitative.l) < length(trait_sum)){
    var.quantitative.l <- c(var.quantitative.l, rep(test,length.out=length(trait_sum)-length(var.quantitative.l)))
  }

  ntraits <- length(trait_sum)
  n.additive <- c(n.additive, rep(0, length.out=ntraits-length(n.additive)))
  n.dominant <- c(n.dominant, rep(0, length.out=ntraits-length(n.dominant)))
  n.equal.additive <- c(n.equal.additive, rep(0, length.out=ntraits-length(n.equal.additive)))
  n.equal.dominant <- c(n.equal.dominant, rep(0, length.out=ntraits-length(n.equal.dominant)))
  n.qualitative <- c(n.qualitative, rep(0, length.out=ntraits-length(n.qualitative)))
  n.quantitative <- c(n.quantitative, rep(0, length.out=ntraits-length(n.quantitative)))

  if(length(unlist(c(var.qualitative.l, var.quantitative.l, var.additive.l, var.dominant.l)))>0){
    ntraits <- max(length(trait_sum), length(var.additive.l),length(var.dominant.l), length(var.qualitative.l), length(var.quantitative.l) )
    n.additive <- c(n.additive, rep(0, length.out=ntraits-length(n.additive)))
    n.dominant <- c(n.dominant, rep(0, length.out=ntraits-length(n.dominant)))
    n.equal.additive <- c(n.equal.additive, rep(0, length.out=ntraits-length(n.equal.additive)))
    n.equal.dominant <- c(n.equal.dominant, rep(0, length.out=ntraits-length(n.equal.dominant)))
    n.qualitative <- c(n.qualitative, rep(0, length.out=ntraits-length(n.qualitative)))
    n.quantitative <- c(n.quantitative, rep(0, length.out=ntraits-length(n.quantitative)))
    trait_sum <- n.additive + n.dominant + n.qualitative + n.quantitative + n.equal.additive + n.equal.dominant
    if(length(var.additive.l) < length(trait_sum)){
      var.additive.l <- rep(var.additive.l, length.out=length(trait_sum))
    }
    if(length(var.dominant.l) < length(trait_sum)){
      var.dominant.l <- rep(var.dominant.l, length.out=length(trait_sum))
    }
    if(length(var.qualitative.l) < length(trait_sum)){
      var.qualitative.l <- rep(var.qualitative.l, length.out=length(trait_sum))
    }
    if(length(var.quantitative.l) < length(trait_sum)){
      var.quantitative.l <- rep(var.quantitative.l, length.out=length(trait_sum))
    }
  }


  if(length(real.bv.dice)>0){
    if(is.list(real.bv.dice)){
      mdepth <- 0
      for(index in 1:length(real.bv.dice)){
        if(is.data.frame(real.bv.dice[[index]][[1]]) || is.matrix(real.bv.dice[[index]][[1]])){
          mdepth <- 1
        }
      }
      if(mdepth==0){
        for(index in 1:length(real.bv.dice)){
          if(length(real.bv.dice[[index]])>0){
            for(index2 in 1:length(real.bv.dice[[index]])){
              if(is.data.frame(real.bv.dice[[index]][[1]][[index2]]) || is.matrix(real.bv.dice[[index]][[1]][[index2]]) ){
                mdepth <- 2
              }
            }
          }
        }
      }
    }
    if(mdepth == 1){
      real.bv.dice <- list(real.bv.dice)
    }
    if(mdepth == 0){
      stop("Illegal input for real.bv.dice")
    }
    for(index in 1:length(real.bv.dice)){
      if(length(real.bv.dice[[index]])>0){
        for(index2 in 1:length(real.bv.dice[[index]][[1]])){
          if(length(real.bv.dice[[index]][[2]][[index2]]) != (3 ^nrow(real.bv.dice[[index]][[1]][[index2]]))){
            warning("")
            stop("Length of effects does not match with involved effect SNPs - should be 3^(effect SNPs) (0..0, 0..01, ..., 2..2)")
          }
        }
      }
    }
  }


  if(length(population)>0){
    if(length(real.bv.add)==0 && replace.traits==FALSE){
      real.bv.add <- population$info$real.bv.add
      if(length(population$info$real.bv.add)>0){
        real.bv.add[[length(population$info$real.bv.add)]] <- NULL
      }

    } else if(replace.traits==FALSE){
      if(!is.list(real.bv.add)){
        real.bv.add <- list(real.bv.add)
      }
      real.bv.add <- c(population$info$real.bv.add, real.bv.add)

      if(length(population$info$real.bv.add)>0){
        real.bv.add[[length(population$info$real.bv.add)]] <- NULL
      }
    }
    if(length(real.bv.mult)==0 && replace.traits==FALSE){
      real.bv.mult <- population$info$real.bv.mult
      if(length(population$info$real.bv.mult)>0){
        real.bv.mult[[length(population$info$real.bv.mult)]] <- NULL
      }

    } else if(replace.traits==FALSE){
      if(!is.list(real.bv.mult)){
        real.bv.mult <- list(real.bv.mult)
      }
      real.bv.mult <- c(population$info$real.bv.mult, real.bv.mult)
      if(length(population$info$real.bv.mult)>0){
        real.bv.mult[[length(population$info$real.bv.mult)]] <- NULL
      }
    }
    if(length(real.bv.dice)==0 && replace.traits==FALSE){
      real.bv.dice <- population$info$real.bv.dice
      if(length(population$info$real.bv.dice)>0){
        real.bv.dice[[length(population$info$real.bv.dice)]] <- NULL
      }

    } else if(replace.traits==FALSE){
      if(!is.list(real.bv.dice)){
        real.bv.dice <- list(real.bv.dice)
      }
      real.bv.dice <- c(population$info$real.bv.dice, real.bv.dice)
      if(length(population$info$real.bv.dice)>0){
        real.bv.dice[[length(population$info$real.bv.dice)]] <- NULL
      }
    }

  }
  if(length(real.bv.add)>0 && !is.list(real.bv.add)){
    real.bv.add <- list(real.bv.add)
  }
  if(length(real.bv.mult)>0 && !is.list(real.bv.mult)){
    real.bv.mult <- list(real.bv.mult)
  }
  if(length(real.bv.dice)>0 && !is.list(real.bv.dice)){
    real.bv.dice <- list(real.bv.dice)
  }

  {
  # Check for missingness in real.bvs and replace with reasonable inputs
  cum_snp <- cumsum(population$info$snp)

  if(length(real.bv.add)>0){
    for(index in 1:length(real.bv.add)){
      while(sum(is.na(real.bv.add[[index]][,1:2]))>0){

        effect_marker <- (1:sum(population$info$snp))
        if(length(exclude.snps)>0){
          effect_marker <- effect_marker[-exclude.snps]
        }


        add_marker <- sample(effect_marker, nrow(real.bv.add[[index]]), replace=if(nrow(real.bv.add[[index]])>length(effect_marker)){TRUE} else{FALSE})
        add_snp <- real.bv.add[[index]][,1]
        add_chromo <- real.bv.add[[index]][,2]

        for(index2 in (1:nrow(real.bv.add[[index]]))[is.na(add_snp) | is.na(add_chromo)]){
          add_chromo[index2] <- sum(add_marker[index2] > cum_snp) + 1
          add_snp[index2] <- add_marker[index2] - c(0,cum_snp)[add_chromo[index2]]
        }

        enter <- add_chromo==real.bv.add[[index]][,2] | is.na(real.bv.add[[index]][,2])

        real.bv.add[[index]][enter,1:2] <- cbind(add_snp, add_chromo)[enter,]
      }
    }
  }

  if(length(real.bv.mult)>0){
    for(index in 1:length(real.bv.mult)){
      for(columns in c(0,2)){
        while(sum(is.na(real.bv.mult[[index]][,1:2+columns]))>0){

          effect_marker <- (1:sum(population$info$snp))
          if(length(exclude.snps)>0){
            effect_marker <- effect_marker[-exclude.snps]
          }


          add_marker <- sample(effect_marker, nrow(real.bv.mult[[index]]), replace=if(nrow(real.bv.mult[[index]])>length(effect_marker)){TRUE} else{FALSE})
          add_snp <- real.bv.mult[[index]][,1]
          add_chromo <- real.bv.mult[[index]][,2]

          for(index2 in (1:nrow(real.bv.mult[[index]]))[is.na(add_snp) | is.na(add_chromo)]){
            add_chromo[index2] <- sum(add_marker[index2] > cum_snp) + 1
            add_snp[index2] <- add_marker[index2] - c(0,cum_snp)[add_chromo[index2]]
          }

          enter <- add_chromo==real.bv.mult[[index]][,2+columns] | is.na(real.bv.mult[[index]][,2])

          real.bv.mult[[index]][enter,1:2+columns] <- cbind(add_snp, add_chromo)[enter,]
        }
      }

    }
  }
  }


  so_far <- max(length(real.bv.dice), length(real.bv.add), length(real.bv.mult))
  if(length(trait_sum)){
    for(index_trait in 1:length(trait_sum)){
      var_additive <- var.additive.l[[index_trait]]
      var_dominant <- var.dominant.l[[index_trait]]
      var_qualitative <- var.qualitative.l[[index_trait]]
      var_quantitative <- var.quantitative.l[[index_trait]]
      if(n.additive[index_trait]>0 && length(var_additive)<n.additive[index_trait]){
        if(length(var_additive)==0){
          var_additive <- 1
        }
        var_additive <- rep(1, length.out=n.additive[index_trait])
        var.additive.l[[index_trait]] <- var_additive
      }
      if(n.dominant[index_trait]>0 && length(var_dominant)<n.dominant[index_trait]){
        if(length(var_dominant)==0){
          var_dominant <- 1
        }
        var_dominant <- rep(1, length.out=n.dominant[index_trait])
        var.dominant.l[[index_trait]] <- var_dominant
      }
      if(n.qualitative[index_trait]>0 && length(var_qualitative)<n.qualitative[index_trait]){
        if(length(var_qualitative)==0){
          var_qualitative <- 1
        }
        var_qualitative <- rep(1, length.out=n.qualitative[index_trait])
        var.qualitative.l[[index_trait]] <- var_qualitative
      }
      if(n.quantitative[index_trait]>0 && length(var_quantitative)<n.quantitative[index_trait]){
        if(length(var_quantitative)==0){
          var_quantitative <- 1
        }
        var_quantitative <- rep(1, length.out=n.quantitative[index_trait])
        var.quantitative.l[[index_trait]] <- var_quantitative
      }

      if(length(var_additive)!= n.additive[index_trait]){
        n.additive[index_trait] <- length(var_additive)
      }
      if(length(var_dominant)!= n.dominant[index_trait]){
        n.dominant[index_trait] <- length(var_dominant)
      }
      if(length(var_qualitative)!= n.qualitative[index_trait]){
        n.qualitative[index_trait] <- length(var_qualitative)
      }
      if(length(var_quantitative)!= n.quantitative[index_trait]){
        n.quantitative[index_trait] <- length(var_quantitative)
      }



      snpdata <- population$info$snp

      #This part is only needed in creating.diploid
      #if(nsnp>0){
      #  snpdata <- c(snpdata, nsnp)
      #} else if(is.matrix(dataset) && nrow(dataset)){
      #  snpdata <- c(snpdata, nrow(dataset))
      #}
      #

      # Generating additive
      effect_marker <- (1:sum(snpdata))
      if(length(exclude.snps)>0){
        effect_marker <- effect_marker[-exclude.snps]
      }

      add_marker <- sample(effect_marker, n.additive[index_trait], replace=if(n.additive[index_trait]>length(effect_marker)){TRUE} else{FALSE})
      dom_marker <- sample(effect_marker, n.dominant[index_trait], replace=if(n.dominant[index_trait]>length(effect_marker)){TRUE} else{FALSE})
      add_marker1 <- sample(effect_marker, n.equal.additive[index_trait], replace=if(n.equal.additive[index_trait]>length(effect_marker)){TRUE} else{FALSE})
      dom_marker1 <- sample(effect_marker, n.equal.dominant[index_trait], replace=if(n.equal.dominant[index_trait]>length(effect_marker)){TRUE} else{FALSE})
      epi1_marker <- sample(effect_marker, n.quantitative[index_trait]*2, replace=if(n.quantitative[index_trait]*2>length(effect_marker)){TRUE} else{FALSE})
      epi2_marker <- sample(effect_marker, n.qualitative[index_trait]*2, replace=if(n.qualitative[index_trait]*2>length(effect_marker)){TRUE} else{FALSE})



      cum_snp <- cumsum(snpdata)
      real.bv.add.new <- NULL
      real.bv.mult.new <- NULL
      if(n.additive[index_trait]>0){
        add_snp <- add_chromo <- numeric(n.additive[index_trait])
        for(index in 1:n.additive[index_trait]){
          add_chromo[index] <- sum(add_marker[index] > cum_snp) + 1
          add_snp[index] <- add_marker[index] - c(0,cum_snp)[add_chromo[index]]
        }
        add_effect <- stats::rnorm(n.additive[index_trait], 0, var_additive)
        real.bv.add.new <- cbind(add_snp, add_chromo, add_effect,0,-add_effect)
      }

      if(n.equal.additive[index_trait]>0){
        add_snp1 <- add_chromo1 <- numeric(n.equal.additive[index_trait])
        for(index in 1:n.equal.additive[index_trait]){
          add_chromo1[index] <- sum(add_marker1[index] > cum_snp) + 1
          add_snp1[index] <- add_marker1[index] - c(0,cum_snp)[add_chromo1[index]]
        }
        add_effect1 <- effect.size.equal.add
        real.bv.add.new <- rbind(real.bv.add.new, cbind(add_snp1, add_chromo1,  -add_effect1, 0, add_effect1))

      }

      if(n.dominant[index_trait]>0){
        dom_snp <- dom_chromo <- numeric(n.dominant[index_trait])
        for(index in 1:n.dominant[index_trait]){
          dom_chromo[index] <- sum(dom_marker[index] > cum_snp) + 1
          dom_snp[index] <- dom_marker[index] - c(0,cum_snp)[dom_chromo[index]]
        }

        dom_effect <- stats::rnorm(n.dominant[index_trait], 0, var_dominant)

        if(dominant.only.positive){
          temp1 <- dom_effect
          temp1[temp1<0] <- 0
        } else{
          temp1 <- dom_effect
        }
        real.bv.add.new <- rbind(real.bv.add.new, cbind(dom_snp, dom_chromo, 0 ,temp1,dom_effect))

      }

      if(n.equal.dominant[index_trait]>0){
        dom_snp1 <- dom_chromo1 <- numeric(n.equal.dominant[index_trait])
        for(index in 1:n.equal.dominant[index_trait]){
          dom_chromo1[index] <- sum(dom_marker1[index] > cum_snp) + 1
          dom_snp1[index] <- dom_marker1[index] - c(0,cum_snp)[dom_chromo1[index]]
        }
        dom_effect1 <- effect.size.equal.dom
        real.bv.add.new <- rbind(real.bv.add.new, cbind(dom_snp1, dom_chromo1, 0 ,dom_effect1, dom_effect1))

      }

      if(n.quantitative[index_trait]){
        epi1_snp <- epi1_chromo <- numeric(n.quantitative[index_trait]*2)
        for(index in 1:(n.quantitative[index_trait]*2)){
          epi1_chromo[index] <- sum(epi1_marker[index] > cum_snp) + 1
          epi1_snp[index] <- epi1_marker[index] - c(0,cum_snp)[epi1_chromo[index]]
        }

        effect_matrix <- matrix(0,nrow=n.quantitative[index_trait], ncol=9)
        for(index in 1:n.quantitative[index_trait]){
          d1 <- sort(abs(stats::rnorm(3, 0, var_quantitative[index])))
          d2 <- sort(abs(stats::rnorm(3, 0, var_quantitative[index])))
          effect_matrix[index,] <- c(d1*d2[1], d1*d2[2], d1*d2[3])
        }
        real.bv.mult.new <- cbind(epi1_snp[1:n.quantitative[index_trait]], epi1_chromo[1:n.quantitative[index_trait]],
                                  epi1_snp[-(1:n.quantitative[index_trait])], epi1_chromo[-(1:n.quantitative[index_trait])],
                                  effect_matrix)
      }

      if(n.qualitative[index_trait]>0){
        epi2_snp <- epi2_chromo <- numeric(n.qualitative[index_trait]*2)
        for(index in 1:(n.qualitative[index_trait]*2)){
          epi2_chromo[index] <- sum(epi2_marker[index] > cum_snp) + 1
          epi2_snp[index] <- epi2_marker[index] - c(0,cum_snp)[epi2_chromo[index]]
        }

        effect_matrix <- matrix(0,nrow=n.qualitative[index_trait], ncol=9)
        for(index in 1:n.qualitative[index_trait]){

          d1 <- -abs(stats::rnorm(9, 0, var_qualitative[index]))
          d1[c(3,7)] <- -d1[c(3,7)]
          effect_matrix[index,] <- d1
        }
        real.bv.mult.new <- rbind(real.bv.mult.new, cbind(epi2_snp[1:n.qualitative[index_trait]], epi2_chromo[1:n.qualitative[index_trait]],
                                                          epi2_snp[-(1:n.qualitative[index_trait])], epi2_chromo[-(1:n.qualitative[index_trait])],
                                                          effect_matrix))
      }

      real.bv.add[[index_trait+so_far]] <- real.bv.add.new
      real.bv.mult[[index_trait+so_far]] <- real.bv.mult.new

    }
  }




  if(length(real.bv.add)>0 && !is.list(real.bv.add)){
    real.bv.add <- list(real.bv.add)
  }
  if(length(real.bv.mult)>0 && !is.list(real.bv.mult)){
    real.bv.mult <- list(real.bv.mult)
  }
  if(length(real.bv.dice)>0 && !is.list(real.bv.dice)){
    real.bv.dice <- list(real.bv.dice)
  }

  nbv <- max(length(real.bv.add), length(real.bv.mult), length(real.bv.dice), if(length(trait_sum)>1){length(trait_sum)} else{0})
  if(nbv >= bv.total){
    bv.total <- nbv
    bv.calc <- nbv
    bv.random <- rep(FALSE, bv.total)
    bv.random.variance <- c(rep(0, nbv))
  }
  if(bv.total > nbv){
    if(length(polygenic.variance)< (bv.total - nbv)){
      polygenic.variance <- rep(polygenic.variance, bv.total - nbv)
    }
    bv.random <- c(rep(FALSE, nbv), rep(TRUE, bv.total - nbv))

    bv.random.variance <- c(rep(0, nbv), polygenic.variance)
    bv.calc <- nbv +1
  }

  population$info$bve <- FALSE
  population$info$bv.calculated <- FALSE
  population$info$breeding.totals <- list()
  population$info$bve.data <- list()
  population$info$bv.nr <- 1 # default um fallunterscheidung zu vermeiden
  population$info$bv.random <- bv.random
  population$info$bv.random.variance <- bv.random.variance

  population$info$phenotypic.transform <- rep(FALSE, bv.total)
  population$info$phenotypic.transform.function <- list()

  store1 <- population$info$is.maternal
  store2 <- population$info$is.paternal
  store3 <- population$info$is.combi
  store4 <- population$info$phenotypic.transform
  store5 <- population$info$phenotypic.transform.function
  store6 <- population$info$bv.random
  store7 <- population$info$bv.random.variance
  store8 <- population$info$bve.mult.factor
  store9 <- population$info$bve.poly.factor
  store10 <- population$info$base.bv

  if(length(is.maternal)==0){
    population$info$is.maternal <- rep(FALSE, bv.total)
  } else{
    if(length(is.maternal)==bv.total){
      population$info$is.maternal <- is.maternal
    } else {
      population$info$is.maternal <- c(population$info$is.materal, rep(is.maternal, length.out = bv.total - length(population$info$is.materal)))
    }

  }
  if(length(is.paternal)==0){
    population$info$is.paternal <- rep(FALSE, bv.total)
  } else{
    if(length(is.paternal)==bv.total){
      population$info$is.paternal <- is.paternal
    } else {
      population$info$is.paternal <- c(population$info$is.paternal, rep(is.paternal, length.out = bv.total - length(population$info$is.paternal)))
    }
  }



  population$info$is.combi <- rep(FALSE, bv.total)


  if(length(bve.mult.factor)==0){
    population$info$bve.mult.factor <- rep(1L, bv.total)
  } else{

    if(length(bve.mult.factor)==bv.total){
      population$info$bve.mult.factor <- bve.mult.factor
    } else {
      population$info$bve.mult.factor <- c(population$info$bve.mult.factor, rep(bve.mult.factor, length.out = bv.total - length(population$info$bve.mult.factor)))
    }
  }

  if(length(bve.poly.factor)==0){
    population$info$bve.poly.factor <- rep(1L, bv.total)
  } else{

    if(length(bve.poly.factor)==bv.total){
      population$info$bve.poly.factor <- bve.poly.factor
    } else {
      population$info$bve.poly.factor <- c(population$info$bve.poly.factor, rep(bve.poly.factor, length.out = bv.total - length(population$info$bve.poly.factor)))
    }

  }
  if(length(base.bv)==0){
    population$info$base.bv <- rep(100L, bv.total)
  } else{

    if(length(base.bv)==bv.total){
      population$info$base.bv <- base.bv
    } else {
      population$info$base.bv <- c(population$info$base.bv, rep(base.bv, length.out = bv.total - length(population$info$base.bv)))
    }
  }


  if(bv.total>0 ||length(real.bv.add)>0  || length(real.bv.mult) >0 || length(real.bv.dice)>0){
    population$info$bve <- TRUE
    if(is.list(real.bv.add)){
      population$info$real.bv.add <- real.bv.add
    } else{
      population$info$real.bv.add <- list(real.bv.add)
    }
    if(is.list(real.bv.mult)){
      population$info$real.bv.mult <- real.bv.mult
    } else{
      population$info$real.bv.mult <- list(real.bv.mult)
    }
    if(is.list(real.bv.dice)){
      population$info$real.bv.dice <- real.bv.dice
    } else{
      if(length(real.bv.dice)>0){
        warning("Invalid input for real.bv.dice!")
      }
      population$info$real.bv.dice <- list(real.bv.dice)
    }

    population$info$bv.nr <- bv.total
    population$info$bv.calc <- bv.calc

    population$info$real.bv.length <- c(length(population$info$real.bv.add),
                                         length(population$info$real.bv.mult),
                                         length(population$info$real.bv.dice))

    population$info$real.bv.add[[nbv+1]] <- "placeholder"
    population$info$real.bv.mult[[nbv+1]] <- "placeholder"
    population$info$real.bv.dice[[nbv+1]] <- "placeholder"



  } else if(preserve.bve){
    population$info$bve <- FALSE
    population$info$bv.nr <- 0
    population$info$bv.calc <- 0
    population$info$real.bv.length <- c(0,0,0)
  }



  if(length(new.residual.correlation)==0 &&
     length(population$info$pheno.correlation)>0 &&
     sum(population$info$pheno.correlation)>sum(diag(population$info$pheno.correlation))){
    if(verbose) cat("Residual correlation has been set to zero since new traits were added ")
  }

  if(length(new.breeding.correlation)==0 &&
     length(population$info$bv.correlation)>0 &&
     sum(abs(population$info$bv.correlation))>sum(diag(population$info$bv.correlation))&&
     sum(population$info$is.combi | !population$info$bv.random) < population$info$bv.nr){
    if(verbose) cat("Genetic correlation between non-QTL traits has been set to zero since new traits were added ")
  }

  store11 <- population$info$pheno.correlation
  store12 <- population$info$bv.correlation
  if(bv.total>0){
    population$info$pheno.correlation <- diag(1L, bv.total)
  }

  if(length(new.residual.correlation)>0){
    population$info$pheno.correlation <- t(chol(new.residual.correlation))
  }
  if(bv.total>0){
    population$info$bv.correlation <- diag(1L, bv.total)
  }
  if(length(new.breeding.correlation)>0){
    population$info$bv.correlation <- new.breeding.correlation
  }


  if(replace.traits==FALSE){
    if(length(store1)>0){
      population$info$is.maternal[1:length(store1)] <- store1
    }
    if(length(store2)>0){
      population$info$is.paternal[1:length(store2)] <- store2
    }
    if(length(store3)>0){
      population$info$is.combi[1:length(store3)] <- store3
    }
    if(length(store4)>0){
      population$info$phenotypic.transform[1:length(store4)] <- store4
    }
    if(length(store5)>0){
      population$info$phenotypic.transform[1:length(store4)] <- store5
    }
    if(length(store6)>0){
      population$info$bv.random[1:length(store6)] <- store6
    }
    if(length(store7)>0){
      population$info$bv.random.variance[1:length(store7)] <- store7
    }
    if(length(store8)>0){
      population$info$bve.mult.factor[1:length(store8)] <- store8
    }
    if(length(store9)>0){
      population$info$bve.poly.factor[1:length(store9)] <- store9
    }
    if(length(store10)>0){
      population$info$base.bv[1:length(store10)] <- store10
    }
    if(length(store11)>0){
      population$info$pheno.correlation[1:nrow(store11), 1:nrow(store11)] <- store11
    }
    if(length(store12)>0){
      population$info$bv.correlation[1:nrow(store12), 1:nrow(store12)] <- store12
    }
  }


  for(generation in 1:nrow(population$info$size)){
    counter <- population$info$size[generation,] + 1
    population$breeding[[generation]][[3]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[1]-1) # estimated breeding value
      population$breeding[[generation]][[4]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[2]-1)
      population$breeding[[generation]][[7]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[1]-1) # real genomic value
      population$breeding[[generation]][[8]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[2]-1)
      population$breeding[[generation]][[9]] <- matrix(NA, nrow= population$info$bv.nr, ncol=counter[1]-1) # phenotype
      population$breeding[[generation]][[10]] <- matrix(NA, nrow= population$info$bv.nr, ncol=counter[2]-1)
      population$breeding[[generation]][[19]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[1]-1) # Reliabilities
      population$breeding[[generation]][[20]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[2]-1)
      population$breeding[[generation]][[21]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[1]-1) # Last applied selection index
      population$breeding[[generation]][[22]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[2]-1)
      population$breeding[[generation]][[27]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[1]-1) # offspring phenotype
      population$breeding[[generation]][[28]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[2]-1)
      population$breeding[[generation]][[29]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[1]-1) # number of offspring used
      population$breeding[[generation]][[30]] <- matrix(0, nrow= population$info$bv.nr, ncol=counter[2]-1)
  }

  if(length(shuffle.traits)==0){
    if(length(shuffle.cor)>0){

      if(ncol(shuffle.cor)==population$info$bv.calc){
        shuffle.traits <- 1:population$info$bv.calc
      } else{
        shuffle.traits <- 1:ncol(shuffle.cor)
        warning(paste0("shuffle.traits not specified! use the first ", ncol(shuffle.cor), " traits"))
      }
    }
  }

  if(length(shuffle.traits)>0){
    if(length(shuffle.traits)==1){
      shuffle.traits <- which(population$info$bv.random==FALSE)
    }

    population <- breeding.diploid(population, verbose = FALSE)
    bvs <- get.bv(population, gen=1)
    scalings <- sqrt(diag(stats::var(t(bvs))))
    for(bvnr in shuffle.traits){
      if(length(population$info$real.bv.add[[bvnr]])>0){
        population$info$real.bv.add[[bvnr]][,3:5] <- population$info$real.bv.add[[bvnr]][,3:5] / scalings[bvnr] * scalings[1]
      }


      if(length(population$info$real.bv.mult[[bvnr]])>0){
        population$info$real.bv.mult[[bvnr]][,5:13] <- population$info$real.bv.mult[[bvnr]][,5:13] / scalings[bvnr] * scalings[1]
      }

      if(length(population$info$real.bv.dice[[bvnr]])>0){
        population$info$real.bv.dice[[bvnr]][[2]] <- population$info$real.bv.dice[[bvnr]][[2]] / scalings[bvnr] * scalings[1]
      }


    }
    population$info$bv.calculated <- FALSE

    LT <- chol(shuffle.cor)
    if(nrow(LT)!=length(shuffle.traits)){
      stop("Dimension of shuffle correlation matrix doesnt work with traits to shuffle")
    } else{

      population$info$bv.correlation[shuffle.traits,shuffle.traits] <- t(LT) %*% LT
      if(sum(abs(population$info$bv.correlation[shuffle.traits,shuffle.traits]- shuffle.cor))>0.0001){
        warning("No covariance matrix for genetic correlation given! Values above diagonal used.")
      }

      store.add <- population$info$real.bv.add
      store.mult <- population$info$real.bv.mult
      store.dice <- population$info$real.bv.dice


      col <- 1
      for(index in shuffle.traits){
        new.add <- new.mult <- new.dice1 <- new.dice2 <- NULL
        row <- 1
        for(index2 in shuffle.traits){
          if(length(store.add[[index2]])>0){
            new.add <- rbind(new.add, store.add[[index2]] %*% diag(c(1,1,rep(LT[row,col],3))))
            zeros <- rowSums(abs(new.add[,3:5,drop=FALSE] ))
            new.add <- new.add[zeros>0,,drop=FALSE]
          }
          if(length(store.mult[[index2]])>0){
            new.mult <- rbind(new.mult, store.mult[[index2]] %*% diag(c(1,1,1,1,rep(LT[row,col],9))))
            zeros <- rowSums(abs(new.mult[,5:13,drop=FALSE]))
            new.mult <- new.mult[zeros>0,,drop=FALSE]
          }
          if(length(store.dice[[index2]])>0){
            if(length(store.dice[[index2]][[1]])>0){
              before <- length(new.dice2)
              new.dice1 <- c(new.dice1,store.dice[[index2]][[1]])
              new.dice2 <- c(new.dice2,store.dice[[index2]][[2]])
              for(index3 in (before+1):length(new.dice2)){
                new.dice2[[index3]] <- new.dice2[[index3]] * LT[row,col]
              }
            }

          }
          row <- row +1
        }

        # DONT REMOVE NULL - MORE WORK NEEDED HERE!
        if(length(new.add)==0){

        } else{
          population$info$real.bv.add[[index]] <- new.add
        }
        if(length(new.mult)==0){

        } else{
          population$info$real.bv.mult[[index]] <- new.mult
        }
        if(length(new.add)==0){

        } else{
          population$info$real.bv.dice[[index]] <- list(new.dice1,new.dice2)
        }

        col <- col +1
      }

    }

    for(index in shuffle.traits){
      population$info$real.bv.length[1] <- max(population$info$real.bv.length[1], if(length(population$info$real.bv.add[[index]])>0){index} else{0})
      population$info$real.bv.length[2] <- max(population$info$real.bv.length[2], if(length(population$info$real.bv.mult[[index]])>0){index} else{0})
      population$info$real.bv.length[3] <- max(population$info$real.bv.length[3], if(length(population$info$real.bv.dice[[index]][[1]])>0){index} else{0})
    }



    for(index in 1:population$info$bv.nr){
      if(length(population$info$real.bv.add[[index]])>0){
        t <- population$info$real.bv.add[[index]]
        take <- sort(t[,1]+ cumsum(c(0,population$info$snp))[t[,2]], index.return=TRUE)
        t <- t[take$ix,,drop=FALSE]
        take <- sort(t[,1]+ t[,2] * 10^10)
        keep <- c(0,which(diff(take)!=0), length(take))
        if(length(keep) <= (nrow(t)+1)){
          for(index2 in 2:(length(keep))){
            t[keep[index2],3:5] <- colSums(t[(keep[index2-1]+1):keep[index2],3:5, drop=FALSE])
          }
          population$info$real.bv.add[[index]] <- t[keep,]
        }
      }
    }


  }

  if(bv.total){
    if(length(population$info$trait.name)>0 & replace.traits==FALSE){
      trait.name <- c(population$info$trait.name, trait.name)
    }
    population$info$trait.name <- trait.name
    if(length(trait.name)<bv.total){
      population$info$trait.name <- c(population$info$trait.name, paste0("Trait ", (length(trait.name)+1):bv.total))
    }
  }

  # Add traits with no generated phenotypes
  temp1 <- rep(0, population$info$bv.nr)
  for(gen in 1:length(population$breeding)){
    for(sex in 1:2){
      if(length(population$breeding[[gen]][[sex]])>0){
        for(index in 1:length(population$breeding[[gen]][[sex]])){
          population$breeding[[gen]][[sex]][[index]][[15]] <- temp1
        }
      }
    }
  }


  if(remove.invalid.qtl && length(population$info$real.bv.add)>1){
    for(index in 1:(length(population$info$real.bv.add)-1)){
      removes <- which(population$info$real.bv.add[[index]][,1] > population$info$snp[population$info$real.bv.add[[index]][,2]])
      if(length(removes)>0){
        population$info$real.bv.add[[index]] <- population$info$real.bv.add[[index]][-removes,,drop=FALSE]
        if(verbose) cat(paste0(removes, " QTL-effects entered on markers that do not exist for ", population$info$trait.name[index], ".\n"))
        if(verbose) cat(paste0(nrow(population$info$real.bv.add[[index]]), " QTL-effects remain.\n"))
      }
    }
    for(index in 1:(length(population$info$real.bv.mult)-1)){
      removes <- which(population$info$real.bv.mult[[index]][,1] > population$info$snp[population$info$real.bv.mult[[index]][,2]])
      if(length(removes)>0){
        population$info$real.bv.mult[[index]] <- population$info$real.bv.mult[[index]][-removes,,drop=FALSE]
        if(verbose) cat(paste0(removes, " QTL-effects entered on markers that do not exist for ", population$info$trait.name[index], ".\n"))
        if(verbose) cat(paste0(nrow(population$info$real.bv.mult[[index]]), " QTL-effects remain.\n"))
      }
    }
    for(index in 1:(length(population$info$real.bv.mult)-1)){
      removes <- which(population$info$real.bv.mult[[index]][,3] > population$info$snp[population$info$real.bv.mult[[index]][,4]])
      if(length(removes)>0){
        population$info$real.bv.mult[[index]] <- population$info$real.bv.mult[[index]][-removes,,drop=FALSE]
        if(verbose) cat(paste0(removes, " QTL-effects entered on markers that do not exist for ", population$info$trait.name[index], ".\n"))
        if(verbose) cat(paste0(nrow(population$info$real.bv.mult[[index]]), " QTL-effects remain.\n"))
      }
    }

  }

  population$info$neff <- list()
  if(length(population$info$real.bv.add)>1){
    for(index in 1:(length(population$info$real.bv.add)-1)){
      if(length(population$info$real.bv.add[[index]])>0){
        population$info$neff[[index]] <- 1:nrow(population$info$real.bv.add[[index]])
      }
    }
  }

  if(bv.standard){
    population <- bv.standardization(population, mean.target = mean.target, var.target = var.target)
  }




  return(population)
}
